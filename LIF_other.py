import numpy as np
import matplotlib.pyplot as plt
import math

class LIF:
    def __init__(self):
        self.C = 1e-6 #(uF)
        self.R = 10e6 #(kOhms?)
        self.e_t = -65
        self.uR = -45
        self.thrs = -50
        self.maxV = -50
        self.v = -65

    def step(self, I, dt, method=0,fire = 0):

        if self.v >= self.thrs and fire:
            #print(self.R* self.C)
            self.v = self.e_t
        else:
            if method == 1:
                self.solve_rk4(dt, I)
            elif method == 0:
                self.solve_euler(dt, I)

            if self.v >= self.thrs:
                self.v = self.maxV

    def solve_euler(self, I, dt):
        dv = self.fu(self.v, I) * dt
        self.v += dv

    def solve_rk4(self, I, dt):
        dv1 = self.fu(self.v, I) * dt
        dv2 = self.fu(self.v + dv1 * 0.5, I) * dt
        dv3 = self.fu(self.v + dv2 * 0.5, I) * dt
        dv4 = self.fu(self.v + dv3, I) * dt
        dv = (1 / 6) * (dv1 + dv2 * 2 + dv3 * 2 + dv4)
        self.v += dv

    def fu(self,v, I):
        return (self.e_t- v + self.R * I) / (self.R * self.C)



def plot(neuron, time, dt, I, method, fire):
    # build the v and u vector
    steps = math.ceil(time / dt)
    v = np.zeros(steps)
    z = np.zeros(steps)
    Q = np.zeros(steps)
    v[0] = neuron.e_t
    vTime = np.arange(0, time, dt, dtype=None)

    for i in range(steps):
        Q[i] = I[i] * 0.35 * ((np.cos(vTime[i]/3)+np.sin(vTime[i]/5)+np.cos(vTime[i]/7)+np.sin(vTime[i]/11)+np.cos(vTime[i]/13))**2)
        neuron.step(dt, Q[i], method, fire)
        v[i] = neuron.v
        #print(v[i],'v[{}]'.format(i))


    plt.plot(vTime, v, color='r', label="current")
    #plt.plot(vTime, Q, color='r', label="current")
    plt.title("Integrate and Fire, Ejercicio 4B")
    plt.xlabel("Tiempo [ms]")
    plt.xlim(xmin=0)
    plt.ylim(ymin=-30)
    plt.ylim(ymax=-70)

    plt.ylabel("Potencial de Membrana")
    plt.show()

def test():
    neuron = LIF()
    time = 200
    dt = 0.5
    steps = math.ceil(time / dt) #10000'
    I = [20e-7 for i in range(steps)]
    print(10e-7)
    plot(neuron, time, dt, I, 1, 1)


test()
