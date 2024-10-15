import sys
import numpy as np

_, filename  = sys.argv

def mean_time(filename):
    with open(filename,'r') as time:
        lines = time.readlines()
        time_list = [x.strip().split() for x in lines]
        res = []
        for j in time_list:
            res.append(float(j[0]))
        print(res)
        print("The average time of the function "+str(len(res)) +" iterations is %.6f sec" % np.mean(res))
    return np.mean(res)

x = mean_time(filename)

