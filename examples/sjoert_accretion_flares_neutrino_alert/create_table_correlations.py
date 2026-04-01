import numpy as np
import os
import pandas as pd
import pickle

analysis_cache_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "cache/"
)

if __name__ == "__main__":
    
    
    with open(os.path.join(
        analysis_cache_dir,"correlations.pkl"
    ), "rb") as fp:
         cor = pickle.load(fp)
            
    print(
        '| i | Run + event number |   Alert   | Accretion flare correlated |  TS   | Energy deposition |'
    )
    
    for i in cor:
        nu = os.path.basename(i[0])
        ic_id = f"IC{nu.split('IceCube-')[-1].split('.')[0]}"
        i.append(ic_id)
        
    cor.sort(key=lambda x: x[2], reverse=True)
    names = [i[3] for i in cor]
    ts = np.array(cor)[:,2]
    evttypes = np.array(cor)[:,3]
    led_ts = np.sum([float(x) for x in ts[evttypes=="LED"]])
    hed_ts = np.sum([float(x) for x in ts[evttypes=="HED"]])
    
    [print(
        f"| {n} | {i[0]} | {names[n]} | {i[1]} | {i[2]:.3e} | {i[3]} |"
    ) for n, i in enumerate(cor)]
    print('-------------------------------------------------------------')
    print(
        f"Test statistic due to LED events: {led_ts:.2f}\n"
        f"Test statistic due to HED events: {hed_ts:.2f}"
    )
