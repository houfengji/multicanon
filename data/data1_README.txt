data1 was made with one planet with following parameters:
K1     = 100 m/s
omega1 = 2pi/100 rad/d,  period = 100 days
phi1   = 1 rad
e      = 0.05
varpi1 = 4 rad

An observation time starts from day 1000, the observation interval dt follows
dt     ~ 1/lambda*exp(-dt/lambda), lambda = 10 days

An observation error, Gaussian i.i.d. with mean 0 and stdev 5 m/s, was added to 
the data, the absolute value of which was recorded in the 3rd column of the
data.

A jitter, Gaussian i.i.d. with mean 0 and stdev 10 m/s, was added to the data.

An offset of 2 m/s was added to the data.
