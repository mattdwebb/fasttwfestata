set type double
insheet using tests/fixtures/data.csv, comma names clear
fasttwfe y x, time(time) cluster(state) gvar(g) inference(native jackknife)
