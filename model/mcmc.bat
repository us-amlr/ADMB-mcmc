krill -ind krill.dat -nox -mcmc 5000000 -mcsave 1000 -nosdmcmc >> out1.txt # runs faster with output to file
krill -ind krill.dat -mceval >> out2.txt                                   # instead of screen

DEL out1.txt /S
DEL out2.txt /S
