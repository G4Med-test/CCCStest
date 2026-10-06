#!/usr/bin/env python

# CCCStest portal parser: metadata() reads the macro, parse() the captured stdout.
import os
import math
import copy
from collections import defaultdict

from geantval import energy_mev, getJSON, one_command, single_run

FRAGMENT_NAMES = [
    "Total",
    "H", "He", "Li", "Be", "B", "C", "N", "O", "alpha", "deuteron", "e-", "et", "gamma", "kaon+", "kaon0L", "kaon0S",
    "lambda", "neutron", "pi", "anti_proton", "proton", "triton"]

ISOTOPES = {'B': [7, 8, 9, 10, 11, 12, 13], 'Be': [5, 6, 7, 8, 9, 10, 11, 12], 'C': [8, 9, 10, 11, 12, 13], 'H': [4, 5, 6, 7], 'He': [3, 5, 6, 7, 8, 9, 10],
            'Li': [4, 5, 6, 7, 8, 9, 10, 11, 13], 'N': [10, 11, 12, 13]}


def read_stdout(filename):
    fragXS = defaultdict(float)
    fragEntries = defaultdict(int)

    fragXS_T = defaultdict(float)
    seen = set()
    fragEntries_T = defaultdict(int)

    with open(filename) as myfile:
        for line in myfile:
            line = line.strip()
            line_l = line.split()
            if len(line_l) > 7 and line_l[3].startswith("MeanFreePath=") and line_l[6].startswith("CrossSection="):
                el = line_l[0]
                if el not in FRAGMENT_NAMES[1:]:  # Skip "Total" just in case
                    print("Element", el, "NOT FOUND when reading isotope crossection")
                    continue
                else:
                    # print "Element", el, "FOUND when reading isotope crossection"
                    try:
                        isot = int("".join(x for x in line_l[1] if x.isdigit()))
                    except ValueError:
                        isot = 0

                    if el in ISOTOPES and isot not in ISOTOPES[el]:
                        print("Isotope", line_l[1], "NOT FOUND in isotope dictionary")
                        continue

                    # print "Isotope", line_l[1], "FOUND in isotope dictionary"
                    fragXS[el + str(isot)] = float(line_l[-2].split('=', 1)[-1])
                    fragEntries[el + str(isot)] = int(line_l[2])

            elif len(line_l) == 4 and line_l[-1] == "[millibarn]":
                el = line_l[0]
                if el.strip() not in FRAGMENT_NAMES:
                    print("Element", el, "NOT FOUND when reading full crossection")
                    continue
                else:
                    # print "Element", el, "FOUND when reading full crossection"
                    if el in seen:
                        raise ValueError("Repeated summary: expected one complete serial run")
                    seen.add(el)
                    fragXS_T[el] = float(line_l[2])
                    fragEntries_T[el] = int(line_l[1])

    if "Total" not in fragXS_T:
        raise ValueError("No CCCStest Total cross section found in stdout")
    fragXS_T["Total"] /= 4.
    return fragXS, fragEntries, fragXS_T, fragEntries_T


PLOT_ALL_FRAGMENTS = False  # Legacy option: also plot every fragment type.


def parse(job):
    filename = os.path.join(job["path"], 'test_stdout.txt')
    print("Parsing", filename)
    fragXS, fragEntries, fragXS_T, fragEntries_T = read_stdout(filename)
    FRAGMENT_LABELS = copy.copy(FRAGMENT_NAMES)
    FRAGMENT_LABELS[0] = "Total * 0.25"

    if PLOT_ALL_FRAGMENTS:
        binEdgeLow = []
        binEdgeHigh = []
        binContent = []
        yStatErrorsPlus = []
        yStatErrorsMinus = []
        for frag in FRAGMENT_NAMES:
            binEdgeLow.append(FRAGMENT_NAMES.index(frag) + 0.5)
            binEdgeHigh.append(FRAGMENT_NAMES.index(frag) + 1.5)
            binContent.append(fragXS_T[frag])
            if fragEntries_T[frag] > 0.:
                errorrel = math.sqrt(float(fragEntries_T[frag])) / float(fragEntries_T[frag])
            else:
                errorrel = 0.
            yStatErrorsPlus.append(fragXS_T[frag] * errorrel)
            yStatErrorsMinus.append(fragXS_T[frag] * errorrel)
        rjson = getJSON(job, "histogram",
                        mctool_name="GEANT4",
                        mctool_model=job["PHYSICS_LIST"],
                        testName=job["TEST"],
                        observableName="fragments production cross section",
                        targetName=job["TARGET"],
                        beamParticle=job["PARTICLE"],
                        beamEnergies=[float(job["ENERGY"])],
                        secondaryParticle="None",
                        binContent=binContent,
                        binLabel=FRAGMENT_LABELS,
                        binEdgeLow=binEdgeLow,
                        binEdgeHigh=binEdgeHigh,
                        yStatErrorsPlus=yStatErrorsPlus,
                        yStatErrorsMinus=yStatErrorsMinus,
                        xAxisName="fragment type",
                        yAxisName="cross section, mb",
                        title="fragments production cross section"
                        )
        yield rjson

    binEdgeLow = []
    binEdgeHigh = []
    binContent = []
    yStatErrorsPlus = []
    yStatErrorsMinus = []
    for i, frag in enumerate(["Total", "Li", "Be", "B"]):
        binEdgeLow.append(i + 0.5)
        binEdgeHigh.append(i + 1.5)
        binContent.append(fragXS_T[frag])
        if fragEntries_T[frag] > 0.:
            errorrel = math.sqrt(float(fragEntries_T[frag])) / float(fragEntries_T[frag])
        else:
            errorrel = 0.
        yStatErrorsPlus.append(fragXS_T[frag] * errorrel)
        yStatErrorsMinus.append(fragXS_T[frag] * errorrel)

    rjson = getJSON(job, "histogram",
                    mctool_name="GEANT4",
                    mctool_model=job["PHYSICS_LIST"],
                    testName=job["TEST"],
                    observableName="fragments production cross section",
                    targetName=job["TARGET"],
                    beamParticle=job["PARTICLE"],
                    beamEnergies=[float(job["ENERGY"])],
                    secondaryParticle="None",
                    binContent=binContent,
                    binLabel=["Total * 0.25", "Li", "Be", "B"],
                    binEdgeLow=binEdgeLow,
                    binEdgeHigh=binEdgeHigh,
                    yStatErrorsPlus=yStatErrorsPlus,
                    yStatErrorsMinus=yStatErrorsMinus,
                    xAxisName="fragment type",
                    yAxisName="cross section, mb",
                    title="fragments production cross section"
                    )
    yield rjson

    for iso in ISOTOPES:
        binLabel = []
        binEdgeLow = []
        binEdgeHigh = []
        binContent = []
        yStatErrorsPlus = []
        yStatErrorsMinus = []

        for i, A in enumerate(ISOTOPES[iso]):
            isoname = iso + str(A)
            binLabel.append(isoname)
            binEdgeLow.append(i + 0.5)
            binEdgeHigh.append(i + 1.5)
            binContent.append(fragXS[isoname])
            if fragEntries[isoname] > 0.:
                errorrel = math.sqrt(float(fragEntries[isoname])) / float(fragEntries[isoname])
            else:
                errorrel = 0.
            yStatErrorsPlus.append(fragXS[isoname] * errorrel)
            yStatErrorsMinus.append(fragXS[isoname] * errorrel)

        rjson = getJSON(job, "histogram",
                        mctool_name="GEANT4",
                        mctool_model=job["PHYSICS_LIST"],
                        testName=job["TEST"],
                        observableName="fragments production cross section",
                        targetName=job["TARGET"],
                        beamParticle=job["PARTICLE"],
                        beamEnergies=[float(job["ENERGY"])],
                        secondaryParticle=iso,
                        binContent=binContent,
                        binLabel=binLabel,
                        binEdgeLow=binEdgeLow,
                        binEdgeHigh=binEdgeHigh,
                        yStatErrorsPlus=yStatErrorsPlus,
                        yStatErrorsMinus=yStatErrorsMinus,
                        xAxisName="fragment type",
                        yAxisName="cross section, mb",
                        title="fragments production cross section"
                        )
        yield rjson


def metadata(commands):
    single_run(commands)
    if one_command(commands, "/testhadr/run/printStat") != "true":
        raise ValueError("Enable printStat to export isotope cross sections")
    particle = one_command(commands, "/gun/particle")
    if particle == "ion":
        ion = one_command(commands, "/gun/ion").split()
        if ion[:2] != ["6", "12"]:
            raise ValueError("This validation is configured for C12 ions")
        particle = "C12"
    return {"TEST": "CCCStest",
            "PHYSICS_LIST": one_command(commands, "/physics/addPhysics"),
            "TARGET": one_command(commands, "/testhadr/det/setMat").removeprefix("G4_"),
            "PARTICLE": particle, "ENERGY": energy_mev(one_command(commands, "/gun/energy"))}
