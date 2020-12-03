
import os
import sys
import numpy as np
import joblib
from skimage import filters
from matplotlib import pyplot as plt
from scipy.signal import *


__file__ = "/home/kai/Ribo-seq-Analysis/module_conservation/uniprot_structure_parser.py" # todo

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from module_supplementary.common_functions import *
from module_conservation.functions import *

genome_base = joblib.load("/Users/kemalinecik/Desktop/genome_base_SHORT.joblib")
# genome_base = joblib.load("/home/kai/KEMALINECIK/out/OUTPUT/protein_structure/genome_base.joblib")



def smooth(x, window_len=11, window='hanning'):
    '''
    Smooth the data using a window with requested size.
    Adapted from:
    http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal

    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this:
    return y[(window_len/2-1):-(window_len/2)] instead of just y.
    '''

    if window_len < 3:  return x

    if x.ndim != 1: raise (Exception('smooth only accepts 1 dimension arrays.'))
    if x.size < window_len:  raise (Exception('Input vector needs to be bigger than window size.'))
    win_type = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    if window not in win_type: raise (Exception('Window type is unknown'))

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')

    # saesha modify
    ds = y.shape[0] - x.shape[0]  # difference of shape
    dsb = ds // 2  # [almsot] half of the difference of shape for indexing at the begining
    dse = ds - dsb  # rest of the difference of shape for indexing at the end
    y = y[dsb:-dse]

    return y

for ind, i in enumerate(genome_base):

    arr = genome_base[i]["footprints"][genome_base[i]['mane_transcript_cds_PA']]
    if sum(arr) > 0:


        thr = np.mean(arr[np.nonzero(arr)])

        smooth_arr_thr = np.array([i if i>=thr else 0 for i in arr])


        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 2.5))

        ax.plot(arr, color='gray', alpha=0.5)

        ax.plot(smooth_arr_thr, color="red")
        plt.title(f"{i}     Count: {sum(arr)}     Thr:{thr}")
        ax.figure.tight_layout()
        fig.savefig(f"_delete_{i}.pdf")
        #plt.show()
        plt.close('all')

for ind, i in enumerate(genome_base):

    a = genome_base[i]["footprints"][genome_base[i]['mane_transcript_cds_PA']]
    if sum(arr) > 0:

        arr = smooth(a, window_len=65)

        thr = find_peaks(smooth(arr))
        prom = peak_prominences(arr, thr[0])
        med_prom = np.median(prom[0])
        peaks = thr[0][prom[0]>med_prom]

        smooth_arr_thr = np.array([i if ind in peaks else 0 for ind, i in enumerate(arr)])




        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 2.5))
        ax.plot(a/max(a)*max(arr), color='gray', alpha=0.5)
        ax.plot(arr, color='gray', alpha=0.5)
        ax.plot(smooth_arr_thr, color="red")
        plt.title(f"{i}     Count: {sum(a)}")
        ax.figure.tight_layout()
        fig.savefig(f"_delete_{i}.pdf")
        #plt.show()
        plt.close('all')

