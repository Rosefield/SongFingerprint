import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wav
import scipy.ndimage as ndimage
import datetime
import argparse
from bisect import bisect_left
import utils

def normalize_samples(samples):
    #input is a 16 bit sample, normalize to [-1, 1]
    scale_const = 2**16
    #if there are multiple tracks, just take the first one
    if(len(samples.shape) > 1):
        samples = samples[:, 0]
    norm_samples = [2*(s / (scale_const + 1.0)) -1  for s in samples]

    return norm_samples

#Returns the value of the Gaussian curve at time t, centered at 0
def gaussian(t):

    return np.exp(-np.pi*(t*t))


def print_image(image, delta):
    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.imshow(image, origin='lower', cmap='jet', interpolation='none', aspect='auto')

    plot.set_xlabel("Sample (delta t = {}s)".format(delta))
    plot.set_ylabel("Frequency (hz)")
    plot.set_title("Spectrograph")

    fig.show()

    if(utils.save_files):
        
        fig.savefig("Spectrogram_{}.png".format(datetime.datetime.now().strftime("%s")))

    raw_input("press enter to continue")

def toDb(data):
    return 20*np.log10(data)


'''
https://en.wikipedia.org/wiki/Gabor_transform


'''
def gabor_transform(samples, rate, delta):

    #Constant at which the gaussian is \le 0.00001
    a = 1.9143

    ts = np.arange(-a, a, 2*a / rate)

    #also of size rate
    window = [gaussian(t) for t in ts]
    i_pad = np.zeros(rate)

    seconds = len(samples) / np.float32(rate)

    segment_hop = int(rate * (delta))
    num_segments = np.int32(np.ceil(len(samples) / (rate *delta))) - 1
    samples = np.concatenate((samples, np.zeros(rate)))

    max_freq = 4000

    samples_by_t = np.empty((num_segments, max_freq), dtype=np.float32)

    for i in xrange(num_segments):
        start = segment_hop * i
        segment = samples[start: start + rate]
        windowed = segment * window
        padded = np.append(windowed, i_pad)
        spectrum = np.fft.fft(padded)
        autopower = np.abs(spectrum * np.conj(spectrum))
        samples_by_t[i, :] = autopower[:max_freq]

    #convert amplitudes to decibels
    db = toDb(samples_by_t)

    #transpose so that time is the x axis
    print_image(np.transpose(db), delta)
    
    return samples_by_t


#hash is of the form f1:f2:deltaf1f2 stored in a 32 bit int
#each frequency is given 12 bits, time is 8 bits
def create_hash(freq1, freq2, diff_time):
    freq1 &= (1 << 12) -1
    freq2 &= (1 << 12) -1
    
    diff_time &= (1 << 8) -1
    

    return (freq1 << 20) | (freq2 << 8) | (diff_time)

#returns timestamp in intervals of 10 milliseconds where sample is the number of sample, and delta is distance between samples
def timestamp(sample, delta):

    return int(100 * sample * delta)

def create_hashes(samples, rate):

    delta = .05
    samples = normalize_samples(samples)

    samples_by_t = gabor_transform(samples, rate, delta)

    #get the samples
    num_samples = samples_by_t.shape[0]
    ts = np.arange(num_samples)

    max_frequencies = np.empty((num_samples), dtype=np.int16)

    for i in ts:
        #skip frequencies under 40, since there is the spike at 0hz
        #divide by 2 because there are rate samples, but we can only get rate/2 frequencies
        freq_max = np.argmax(samples_by_t[i, 80:], axis=0) / 2 + 40
        max_frequencies[i] = int(freq_max)

    data = []
    data.append(("Max Frequency", range(num_samples), max_frequencies))
    utils.plot(data, xlabel="Sample", ylabel="Frequency (hz)", title="Max Frequency by Sample", linestyle='')

    fanout = 10

    #Hashes are of the form ( aboslute timestamp in ms, freq1:freq2:deltaf1f2 ) 
    hashes = []

    for i, freq in enumerate(max_frequencies[ : -fanout]):
        in_target_range = lambda x: (x[1] + 500 >= freq or x[1] - 500 <= freq) and x[0] % 3 == 0
        targets = filter(in_target_range, enumerate(max_frequencies[i:]))[:fanout]
        for x, n_freq in targets:
            h = (timestamp(i, delta), 
                    create_hash(freq, n_freq, timestamp(x, delta)))
            hashes.append(h)

    hashes = sorted(hashes, key=lambda x: x[1])

    return hashes

def read_hashes(filename):

    with open(filename) as f:
        hashes = [tuple(map(int, x.split(','))) for x in f]
        return hashes        

#assumes that the sets are already sorted
def matching_hashes(set1, set2):
    
    matches = []
    
    keys = [x[1] for x in set1]
    size = len(set1)

    for needle in set2:
        n_key = needle[1]

        i = bisect_left(keys, n_key)
        if i < size:
            matches.append((set1[i], needle))   
       

    data = []
    x_vals = [x[0][0] for x in matches]
    y_vals = [x[1][0] for x in matches] 
    data.append(("Matches by file times", x_vals, y_vals))

    file = "Constellation_map_{}.png".format(datetime.datetime.now().strftime("%s"))

    utils.plot(data, xlabel="Source File sample", ylabel="Check File sample", linestyle='', filename=file) 


    return matches


def hashes_match(set1, set2):

    matches = matching_hashes(set1, set2)
    
    check_size = len(set2)
    match_size = len(matches)

    #If there are no matches, obviously it isn't the same or most of the hashes don't match
    if(match_size == 0 or float(match_size)/check_size < .5):
        return False


    buckets = {}
    #each bucket is for .2 a second, time stamps are in 10s of milliseconds
    bucket_size = 20.
    round_bucket = lambda x: int(round(x/bucket_size) * bucket_size)

    for match in matches:
        #difference in sample, we want to see if the samples are sequential (i.e. difference is constant)
        diff = match[0][0] - match[1][0]

        rounded = round_bucket(diff)

        if rounded in buckets.iterkeys():
            buckets[rounded] += 1
        else:
            buckets[rounded] = 1

    max_item = max(buckets.iteritems(), key=lambda x: x[1])

    val = max_item[1]
    num_matches = len(matches)

    print "max sequential: {}, total matches: {}, total samples: {}".format(val, num_matches, check_size)
    #Check that the many of the samples are sequential 
    if(float(val)/num_matches > .4):
        print "Found match at time {}s with {} matching samples".format(max_item[0] * .01, val)
        return True, max_item[0]
    else:
        return False, None


def create_and_save_hashes(filename, save=False):
    rate, samples = wav.read(filename)
    h = create_hashes(samples, rate)
    
    if(save):
        name = filename.split(".")[0] + ".hashes"
        with open(name, 'w') as f:
            for x in h:
                f.write("{},{}\n".format(x[0],x[1]))


    return h

def files_match(file1, file2, save= False):

    files = [file1, file2]

    hashes = []

    for file in files:
        if(file.endswith(".wav")):
            h = create_and_save_hashes(file, save)
            hashes.append(h)
        else :
            hashes.append(read_hashes(file))

    isMatch = hashes_match(hashes[0], hashes[1])

    return



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("files", nargs='+')
    parser.add_argument("--save", action="store_true")

    args = parser.parse_args()

    save_hashes = False
    if "save" in args:
        utils.save_files = args.save
        save_hashes = True

    files = args.files
    if len(files) == 1:
        create_and_save_hashes(files[0], save=save_hashes)

    elif len(files) > 1:
        files_match(files[0], files[1], save=save_hashes)

    return


if __name__ == "__main__":
    main()
