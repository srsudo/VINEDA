import math
import obspy
from obspy.core.trace import Trace, Stats
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import scipy.signal
from scipy.signal import butter, lfilter, firwin, filtfilt, medfilt, hilbert

def detect_explosions(x, fs, flow, fhigh, l, nc, beta):

    """
    Algorithm to detect volcanic explosions from a continuous stream of infrasonic data.
    
    :param x: Data samples
    :param fs: Original Sampling rate (Hz)
    :param flow: Lower cutoff frequency of the bandpass filter (Hz)
    :param fhigh: Upper cutoff frequency of the bandpass filter (Hz)
    :param l: Average duration in seconds of the explosions to be detected
    :param nc: Number of sub-bands used in the filter bank
    :param beta: Penalty factor for non-impulsive onset
    :return: (CF, fsp), tuple containing the characteristic function and sampling frequency.
    """

    ################
    # TEST ISAAC ###
    ################
    resultados_isaac = scipy.io.loadmat("Resultados.mat")
    testeo = np.asarray(resultados_isaac["xrd"])

    # Downsampling and decimate the signal

    fsp = math.ceil(2*fhigh)

    # Detrend and decimation
    x.detrend(type='demean')
    data = np.asarray(x.data)

    # F1
    xr = scipy.signal.resample_poly(data, fsp, fs)
    # Adaptive de-trending filter to remove trends in order to avoid spike noise
    n = int(1 * l * fsp)
    w_n = float(flow) / float((fsp / 2.0))

    B = firwin(n, w_n)
    trend = filtfilt(B, 1.0, xr)
    xrd = xr - trend

    # given the issue of parity by the Kernel of scipy, we check the parity of the kernel-
    med_filt_order = int(round((L/10)*fsp))

    if (med_filt_order % 2 ) == 0:
        med_filt_order = med_filt_order -1

    #  Median Filter to avoid spiky signals or instrumental noise
    y = medfilt(xrd, med_filt_order)

    # no edge delay window
    lh = round(l*fsp)
    th = np.arange(0, lh/float(fsp), 1/fsp)
    haux = np.linspace(1, 0, lh)
    k = np.arange(1, int(nc)+1, 1)  # We should keep it one as is freq.

    fl = flow + (k-1)*(fhigh-flow)/float(nc)
    fh = fl + (fhigh-flow) / float(nc)

    v = np.zeros((y.shape[0], int(nc)))
    ei = np.zeros((y.shape[0], int(nc)))

    for m in range(int(nc)):
        fc = (fl[m]+fh[m])/2.0
        h1 = haux * np.sin(2 * np.pi * fc * th)
        v[:, m] = lfilter(h1, 1, y)
        ei[:, m] = np.abs(hilbert(v[:, m]))

    e = np.sum(ei, axis=1)

    # Discriminant detector
    LB = np.round(l*fsp)

    l1 = np.arange(1, 41, 1)
    l2 = [beta * (x - (LB+1)) for x in l1]
    h2 = np.concatenate((l1, l2))
    h2 /= np.linalg.norm(h2)

    # CF function building
    cf = lfilter(h2, 1.0, e, axis=0)

    M = int(math.floor(h2.shape[0]/2.0))  # group delay

    # check this
    cf = np.vstack((cf[M:-1].reshape(-1, 1),
                    np.zeros((M, 1))))


    cf[0:M-1] = 0
    cf[cf <= 0] = 0

    tw = int(round(5*fsp))  # window of 30 seconds to estimate the background noise level

    h3 = (1.0/tw) * np.ones((tw,))

    #filter
    n = lfilter(h3, 1.0, e)  

    M = int(math.floor(tw/2.0))  

    n = np.vstack((n[M:-1].reshape(-1, 1),
                    np.zeros((M, 1))))

    #n[0:M-1] = n[M]
    n[(n.shape[0]-M):] = n[n.shape[0]-M-1]

    #  perform the normalization
    cf[cf != 0] = cf[cf != 0] / (n[cf != 0])

    return (cf, fsp)


def read_mat_file(filename):
    """
    :param filename: String with the filename we want to read.
    :return: the data on the stream
    """
    return scipy.io.loadmat(filename)


def plot_stream(stream, fsp, stream_name):
    """
        Auxiliary function to plot streams of data.

        :param stream: A numpy array containing the stream data.
        :param stream_name: A string containing the stream information
        :return:
        """
    time_vector = np.linspace(0, len(stream)/fsp, num=len(stream))
    plt.title(stream_name)
    plt.plot(time_vector, stream)
    plt.show()
    plt.close()

def convert_to_obspy(signal, fsp):
    """
    Function to convert a numpy array from any infrasound reading into Obspy format.
    This could help us to take advantage of Obspy optimization routines.
    :param signal: A numpy array containg the signal we want to compute.
    :param fsp: The sampling frequency of the signal.
    :return:
    """
    # we substract the mean

    stats = Stats()
    stats.sampling_rate = float(fsp)
    stats.npts = signal.shape[0]
    return Trace(data=signal.reshape(signal.shape[0], ), header=stats)

if __name__=='__main__':

    #########################
    # PARAMETERS DEFINITION #
    #########################

    fs = 100.0    # Original sampling rate
    f_high = 5.0  # Corner top frequency of the filter bank
    f_low = 0.1   # Minimum frequency of the band pass filter
    L = 4.0       # Average duration (s) of the event
    Nc = 3.0      # Number of central frequency used in the filter bank
    beta = 3.0    # Penalty factor for non - impulsive onsets

    # read the data
    x = read_mat_file("PruebaSignal_1.mat")

    # VINEDA is agnostic-free algorithm that can take data from any Python
    # pipeline. In this example, we use a data file directly imported from
    # Matlab, but it can extended to any of the Obspy reading formats.

    signal = np.asarray(x["signal"])
    fs = float(x["fs"])

    #we convert to obspy to take advantage of Obspy processing routines.
    stream = convert_to_obspy(signal, fs)

    results = detect_explosions(stream, stream.stats.sampling_rate, f_low, f_high, L, Nc, beta)

    #results is a tuple containing the cf and the fsp, we can plot it-
    plot_stream(results[0], results[1], "CF")

