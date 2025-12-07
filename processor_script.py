import uproot
import json

import os
import argparse
from coffea.nanoevents import NanoEventsFactory
from topcoffea.scripts.make_html import make_html
from topcoffea.modules import utils
import topcoffea.modules.quad_fit_tools as qft

from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema
import awkward as ak
import numpy as np

import glob
from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
from coffea import util
import vbfprocessor


def main():
    indir = (
        "/eos/uscms/store/user/jennetd/vbf-eft/"
        "VBF_SMEFTsim_topU3l_Direct_cHWtil_HT2/"
        "VBF_SMEFTsim_topU3l_Direct_cHWtil_HT2/241004_221644/0000"
    )

    files = sorted(glob.glob(os.path.join(indir, "*.root")))
    print(f"Found {len(files)} ROOT files in {indir}")
    if len(files) == 0:
        raise RuntimeError("No root files found.")

    sample_name = "VBF_SMEFTsim_topU3l_Direct_cHWtil_HT2"
    fileset = {sample_name: files}

    p = vbfprocessor.VBFProcessor(isMC=True)

    executor = FuturesExecutor(workers=4, status=True)

    runner = Runner(
        executor=executor,
        savemetrics=True,
        schema=NanoAODSchema,
        chunksize=50000,
    )

    output, metrics = runner(
        fileset=fileset,
        processor_instance=p,
        treename="Events",
    )

    os.makedirs("coffea", exist_ok=True)
    outfile = f"coffea/{sample_name}.coffea"
    util.save(output, outfile)

    print(f"Saved {outfile}")
    print("Metrics:", metrics)


if __name__ == "__main__":
    main()
