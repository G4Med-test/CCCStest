# Portal JSON parser

Imported from `geant-config-generator/tests/geant4/CCCStest/parser.py` at commit
`0bf42953c9c98ac6ecccf077e7c7fe2c4ceb7307`, with queue discovery replaced by the
shared ci-workflows exporter.

`config.json` selects `test.mac`. The workflow captures stdout in `test_stdout.txt`;
this is the parser input, not an XML file. The parser derives model, target, beam
and energy from the macro. One serial run is expected per macro. Enable
`/testhadr/run/printStat true` to include isotope cross sections; this is now set
in the validation macro. The unit macro uses the same valid application commands
with one event.

The legacy eight plots, isotope selection, `Total * 0.25` bin and errors
`cross_section / sqrt(entries)` are preserved. Missing production channels remain
zero as in the original parser, but missing or repeated summary tables are errors.

The pipeline uploads `validation-json/results.json` (an array of plot objects)
and `plots/*.json` (the original one-object-per-file format), plus raw artifacts.
See [ci-workflows](https://github.com/G4Med-test/ci-workflows) for local CLI usage.
