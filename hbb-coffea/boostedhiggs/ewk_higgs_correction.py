from __future__ import annotations

from pathlib import Path
import numpy as np
import awkward as ak
import correctionlib

_THIS_DIR = Path(__file__).resolve().parent
_EWK_JSON = _THIS_DIR / "data" / "EWHiggsCorrections.json"


if not _EWK_JSON.exists():
    raise FileNotFoundError(
        f"[EWK] Can't find EW Higgs json at {_EWK_JSON}\n"
        f"Check that json exists or change the relative path."
    )

hew_kfactors = correctionlib.CorrectionSet.from_file(str(_EWK_JSON))

def add_HiggsEW_kFactors(weights, genpart, dataset):
    print("adding EWK")
    """EW Higgs corrections"""

    def get_hpt():
        boson = ak.firsts(
            genpart[
                (genpart.pdgId == 25)
                & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
            ]
        )

        hpt = np.array(ak.fill_none(boson.pt, 0.))
        print("hpt min & max:", float(np.min(hpt)), float(np.max(hpt)))

        return hpt

    if "VBF" in dataset:
        hpt = get_hpt()
        ewkcorr = hew_kfactors["VBF_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        print("correction's min & max:", float(np.min(ewknom)), float(np.max(ewknom)))
        weights.add("VBF_EW", ewknom)

        print("added VBF_EW")

        return hpt
#        hpt = get_hpt()
#        ewkcorr = hew_kfactors["VBF_EW"]
#        ewknom = ewkcorr.evaluate(hpt)
#        weights.add("VBF_EW", ewknom)

    if ("WplusH" in dataset) or ("WminusH" in dataset) or ("ZH" in dataset) or ("VH" in dataset):
        hpt = get_hpt()
        ewkcorr = hew_kfactors["VH_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        weights.add("VH_EW", ewknom)

    if "ttH" in dataset:
        hpt = get_hpt()
        ewkcorr = hew_kfactors["ttH_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        weights.add("ttH_EW", ewknom)
