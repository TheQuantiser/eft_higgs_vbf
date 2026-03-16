import time

import coffea.processor as processor
import hist
from coffea.analysis_tools import PackedSelection, Weights
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
from coffea.nanoevents.methods import nanoaod

NanoAODSchema.warn_missing_crossrefs = False

import pickle
import re

import awkward as ak
import numpy as np
import pandas as pd
import json
import os

#ewk import
from boostedhiggs.ewk_higgs_correction import add_HiggsEW_kFactors

#systematic import
from boostedhiggs.corrections import (
    add_ps_weight,
    add_pdf_weight,
    add_scalevar_7pt,
    add_scalevar_3pt
)

import importlib.util
import boostedhiggs.corrections as corrections

importlib.reload(corrections)

def weight_names(reweight_card):

    with open('dictionary.json') as f:
        d = json.load(f)

    array = []
    weight_index = 0

    with open(reweight_card, "r+") as file1:
        for line in file1.readlines():
            
            if 'launch' in line:
                weight_index += 1
                continue
                
            if 'rwgt' in line:
                continue
                
            array += [[weight_index]+line.split()[1:]]

    df = pd.DataFrame(array,columns=['weight_index','model','param_index','value'])
    df['value'] = df['value'].astype('float')
    df = df[df['value']>0]

    df["model"] = df["model"].str.lower()
    df["param_index"] = df["param_index"].astype(str)
    
    param_name = []
    for index, row in df.iterrows():
        param_name += [d[row['model']][row['param_index']]]

    df['param_name'] = param_name    
    df['weight_name'] = df['param_name']+"="+df['value'].astype('string')

    names = []
    for i in np.unique(df['weight_index']):
        name = ""
        for n in df[df['weight_index']==i]["weight_name"]:
            if len(name)>0:
                name += ","
            name += n
        names += [name]

    return names

def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out

# Look at ProcessorABC to see the expected methods and what they are supposed to do
class VBFProcessor(processor.ProcessorABC):
    def __init__(self, isMC=False, ewkHcorr=False, systematics=False):

        #adding Higgs electroweak correction
        self._ewkHcorr = ewkHcorr
        #adding systematics
        self._systematics = systematics
        
        ################################
        # INITIALIZE COFFEA PROCESSOR
        ################################
            
        ak.behavior.update(nanoaod.behavior)

        q1pt_axis = hist.axis.Regular(200, 0, 1000, name="q1pt", label=r"Quark 1 $p_{T}$ [GeV]")
        q2pt_axis = hist.axis.Regular(200, 0, 1000, name="q2pt", label=r"Quark 2 $p_{T}$ [GeV]")
        hpt_axis = hist.axis.Regular(200, 0, 1000, name="hpt", label=r"Higgs $p_{T}$ [GeV]")
        detaqq_axis = hist.axis.Regular(200, -6, 6, name="detaqq", label=r"$\Delta\eta_{qq}$")
        dphiqq_axis = hist.axis.Regular(200, -3.15, 3.15, name="dphiqq", label=r"$\Delta\phi_{qq}$")
        mqq_axis = hist.axis.Regular(200, 0, 5000, name="mqq", label=r"$m_{qq}$ [GeV]")
        wc_axis = hist.axis.StrCategory([], name="wc", label="WC point", growth=True)
        htxs_axis = hist.axis.IntCategory([], name = "htxs_stage2", label="HTXS stage1.2 code", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label="Systematic", growth=True)

        self.make_output = lambda: { 
            # Test histogram; not needed for final analysis but useful to check things are working
            "q1pt": hist.Hist(
                q1pt_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "q2pt": hist.Hist(
                q2pt_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "hpt": hist.Hist(
                hpt_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "detaqq": hist.Hist(
                detaqq_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "dphiqq": hist.Hist(
                dphiqq_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "mqq": hist.Hist(
                mqq_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            "htxs": hist.Hist(
                htxs_axis,
                wc_axis,
                storage=hist.storage.Weight() #variances are stored in here. can do h.values() and h.variance()
            ),
            "templates": hist.Hist(
                syst_axis,
                htxs_axis,
                wc_axis,
                storage=hist.storage.Weight()
            ),
            
            "EventCount": processor.value_accumulator(int),
            "sumw_all_noEW": processor.value_accumulator(float),
            "sumw_VBF_EW": processor.value_accumulator(float),
            "n_evts": processor.value_accumulator(int),
        }
        
    def process(self, events):
        output = self.make_output()

        ##################
        # OBJECT SELECTION
        ##################

        # Nothing for now, only truth particles

        outgoing = events.LHEPart[events.LHEPart.status==1]
        higgs = ak.firsts(outgoing[outgoing.pdgId == 25])
        quarks = outgoing[outgoing.pdgId<=6]

        q1 = ak.firsts(quarks[:,0:1])
        q2 = ak.firsts(quarks[:,1:2])

        mqq = (q1+q2).mass
        detaqq = q1.eta - q2.eta
        dphiqq = q1.delta_phi(q2)

        htxs_stage2 = events.HTXS["stage1_2_fine_cat_pTjet25GeV"]
        if len(events) > 0:
            htxs_arr = ak.to_numpy(htxs_stage2)
            print("HTXS sample:", htxs_arr[:10])
            print("HTXS unique (first 20):", np.unique(htxs_arr)[:])

        #####################
        # EVENT SELECTION
        #####################

        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        # Nothing for now, everything from generator

        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights = Weights(len(events), storeIndividual=True)
        weights.add('genweight', events.genWeight)

        output["sumw_all_noEW"] += float(ak.sum(events.genWeight))
        output["EventCount"] += int(len(events))
        output["n_evts"] += int(len(events))

       #ew correction diagnostics
        dataset = events.metadata.get("dataset", "")
        
        #weight before EW correction
        w_before = weights.weight()
        
        if self._ewkHcorr:
            try:
                hpt = add_HiggsEW_kFactors(weights, events.GenPart, "VBF")
                print("EWK applied")
            except Exception as e:
                print("EWK exception:", repr(e))

        if self._systematics:
            if "LHEScaleWeight" in events.fields:
                add_scalevar_3pt(weights, events.LHEScaleWeight)
            else:
                add_scalevar_3pt(weights, [])
        print("type:", ak.type(events.LHEScaleWeight))
        print("LHEScaleWeight: ", events.LHEScaleWeight)
        print("variations:", list(weights.variations))

        systematics = [None] + list(weights.variations)

        # SM templates
        for syst in systematics:
            sname = "nominal" if syst is None else syst
            w = weights.weight() if syst is None else weights.weight(modifier=syst)
        
            # Yield-only version:
            output["templates"].fill(
                systematic=sname,
                htxs_stage2=htxs_stage2,
                wc="SM",
                weight=w,
            )
        
        output["sumw_VBF_EW"] += float(ak.sum(weights.weight()))
        
        #weight after EW correction
        w_after = weights.weight()
        
        print("avg weight before:", float(ak.mean(w_before)))
        print("avg weight after: ", float(ak.mean(w_after)))

        k = weights.partial_weight(include=["VBF_EW"])

        is_bin1 = (htxs_stage2 == 221) | (htxs_stage2== 222)
        is_bin2 = (htxs_stage2 == 223) | (htxs_stage2 == 224)
        
        print("avg correction applied to bin1:", float(ak.mean(k[is_bin1])))
        print("avg correction applied to bin2:", float(ak.mean(k[is_bin2])))
        
        w_nom = weights.weight()
        
        y1_nom = float(ak.sum(w_nom[is_bin1]))
        y2_nom = float(ak.sum(w_nom[is_bin2]))

        def safe_ratio(num, den):
            return num / den if den != 0 else float("nan")
            
        for systematic in ["scalevar_3ptUp", "scalevar_3ptDown"]:
            w_systematic = weights.weight(modifier=systematic)
            y1 = float(ak.sum(w_systematic[is_bin1]))
            y2 = float(ak.sum(w_systematic[is_bin2]))
            print(systematic,
              "bin1 ratio:", safe_ratio(y1, y1_nom),
              "bin2 ratio:", safe_ratio(y2, y2_nom))


        ###################
        # FILL HISTOGRAMS
        ###################
        
        try:
            print('try1')
            print(os.getcwd())
            # find the event weight to be used when filling the histograms
            names = weight_names("VBF_SMEFTsim_topU3l_NP1_reweight_card.dat")
            print("names length:", len(names))
            print("first few names:", names[:5])

            # Loop over weights
            for i,n in enumerate(names):

                output['q1pt'].fill(q1pt=q1.pt,
                                    wc=names[i],
                                    weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                   )
                output['q2pt'].fill(q2pt=q2.pt,
                                    wc=names[i],
                                    weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                   )
                output['hpt'].fill(hpt=higgs.pt,
                                   wc=names[i],
                                   weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                  )
                output['detaqq'].fill(detaqq=detaqq,
                                      wc=names[i],
                                      weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                     )
                output['dphiqq'].fill(dphiqq=dphiqq,
                                      wc=names[i],
                                      weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                     )
                output['mqq'].fill(mqq=mqq,
                                   wc=names[i],
                                   weight=weights.weight()*events.LHEReweightingWeight[:,i]
                                  )
                output['htxs'].fill(htxs_stage2=htxs_stage2,
                                   wc=names[i],
                                   weight = weights.weight()*events.LHEReweightingWeight[:,i]
                                   )

            # Also do SM
            output['q1pt'].fill(q1pt=q1.pt,
                                wc='SM',
                                weight=weights.weight()
                                )
            output['q2pt'].fill(q2pt=q2.pt,
                                wc='SM',
                                weight=weights.weight()
                               )
            output['hpt'].fill(hpt=higgs.pt,
                               wc='SM',
                               weight=weights.weight()
                               )
            output['detaqq'].fill(detaqq=detaqq,
                                  wc='SM',
                                  weight=weights.weight()
                                  )
            output['dphiqq'].fill(dphiqq=dphiqq,
                                  wc='SM',
                                  weight=weights.weight()
                                 )
            output['mqq'].fill(mqq=mqq,
                               wc='SM',
                               weight=weights.weight()
                           )
            output['htxs'].fill(htxs_stage2=htxs_stage2,
                                wc='SM',
                                weight=weights.weight(),
                                )
            
            # End if LHE weight
                
        except Exception as e:
            print("except:", repr(e))

            output['q1pt'].fill(q1pt=q1.pt,
                                wc='SM',
                                weight=weights.weight()
                               )
            output['q2pt'].fill(q2pt=q2.pt,
                                wc='SM',
                                weight=weights.weight()
                               )
            output['hpt'].fill(hpt=higgs.pt,
                               wc='SM',
                               weight=weights.weight()
                               )
            output['detaqq'].fill(detaqq=detaqq,
                                  wc='SM',
                                  weight=weights.weight()
                                  )
            output['dphiqq'].fill(dphiqq=dphiqq,
                                  wc='SM',
                                  weight=weights.weight()
                                 )
            output['mqq'].fill(mqq=mqq,
                               wc='SM',
                               weight=weights.weight()
                              )
            output['htxs'].fill(htxs_stage2=htxs_stage2,
                                wc='SM',
                                weight=weights.weight(),
                                )
            
        # End if no LHE weight
    
        output["EventCount"] = len(events)

        return output

    def postprocess(self, accumulator):
        return accumulator
