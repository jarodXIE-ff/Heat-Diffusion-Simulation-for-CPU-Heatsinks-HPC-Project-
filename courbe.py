import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
import numpy as np

_, filename  = sys.argv


def get_all_time(filename):
    with open(filename,'r') as time:
        lines = time.readlines()
        time_list = [x.strip().split() for x in lines]
    print(time_list)
    return time_list

all_time =  get_all_time(filename)


mon_dictionnaire = {}

for valeur, groupe in all_time:

    if groupe in mon_dictionnaire:
        mon_dictionnaire[groupe].append(float(valeur))
    else:
        mon_dictionnaire[groupe] = [float(valeur)]

for cle in mon_dictionnaire.keys():
    mon_dictionnaire[cle] = np.mean(mon_dictionnaire[cle])
print(mon_dictionnaire)

x = []
y = []

for cle, val in mon_dictionnaire.items():
    x.append(cle)
    y.append(val)

plt.plot(x,y, marker='o', linestyle='-')
plt.xlabel("Number of Nodes")
plt.ylabel("Time in sec")
plt.title("Time per nodes for " + filename)
plt.show()
name_courbe = "courbe"+filename+".png"
plt.savefig(name_courbe)