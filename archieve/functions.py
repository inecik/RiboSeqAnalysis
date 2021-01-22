import numpy as np


def smooth(x, window_len=15, window='hanning'):
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

def gauss_kern(size, sizey=None):
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None):
    g = gauss_kern(n, sizey=ny)
