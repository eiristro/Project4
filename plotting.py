from scitools.std import *
import time

T = linspace(0, 1, 101)

ofile1 = open('fEuler.dat', 'r')
ofile2 = open('bEuler.dat', 'r')
ofile3 = open('CN.dat', 'r')

for line1, line2, line3 in zip(ofile1, ofile2, ofile3):
    data1 = line1.split(' ')
    data2 = line2.split(' ')
    data3 = line3.split(' ')

    for i in range(len(data1)-1):
        data1[i] = float(data1[i])
        data2[i] = float(data2[i])
        data3[i] = float(data3[i])
    x = linspace(0, 1, len(data1) - 2)
    t = data1[-2]
    data1 = data1[0:-2]
    data2 = data2[0:-2]
    data3 = data3[0:-2]
    if (t == 0.1) or (t == 0.5):
        figure()
        cla()
        hold('on')
        plot(x, data1, 'r')
        plot(x, data2, 'b')
        plot(x, data3, 'g')

        #Find the analytical solution for comparison
        n = 1
        v = 0
        d = 1
        while d > 0.001:
            dv = (-2/(n*pi))*sin(n*pi*x)*exp(-(n**2)*(pi**2)*t)
            v += dv
            d = sum(dv) / sum(v)
            n += 1
        u = v + (1 - x)
        plot(x, u, 'k')

        legend(['Forward Euler', 'Backward Euler', 'Crank-Nicolson', 'Analytical solution'])
        xlabel('Relative distance from presyntaptic neuron in units of the width of the synaptic cleft')
        ylabel('Neurotransmitter density')
        title('Time: %f' % t)

raw_input('enter to close')
