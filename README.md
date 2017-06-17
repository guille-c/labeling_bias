# Labeling Bias

This repository includes code for calculating how biased data sets are in terms of observational parameters. 

## Usage: 
```
python calculate_bias.py [-h] [--number_objects N] [--int_pars INT_PARS]
                         [--obs_pars OBS_PARS]
                         [--pbb_thresholds PBB_THRESHOLDS [PBB_THRESHOLDS ...]]
                         [--no_zeros] [--bins_obs BINS_OBS]
                         [--log2_bins_int LOG2_BINS_INT] [--N_iter N_ITER]
                         [--labels LABELS [LABELS ...]]
                         table_file label_name [label_name ...]
```

arguments:
- table_file: table with variables and labels.
- label_name: name of the column containing the labels.

- -h, --help: show this help message and exit
- --number_objects: number of objects per bin to be used for calculating the bias.
- --int_pars: Name of the columns used as intrinsic parameters.
- --obs_pars: Name of the columns used as observable parameters.
- --pbb_thresholds: Threshold(s) on the probabilities-
- --no_zeros: Do not consider labels that don't match the pbb. thresholds
- --bins_obs: log2 bins in intrinsic parameters.
- --log2_bins_int: Bins in observable parameters.
- --N_iter: Number of calculations of L to calculate means and standard deviations.
- --labels: Labels to be used.

## Example

If you want to calculate the bias of the catalog at GZ_DR7_all.fit, using 58 objects, a threshold of 0.8 in the probabilities defined by columns p_cs and p_el, 8x8x8 bins in the intrinsic parameters, and 8 bins in the observable parameters you should run it as:

```
python calculate_bias.py GZ_DR7_all.fit p_cs p_el --pbb_thresholds 0.8 0.8 --number_objects 58 --log2_bins_int 9 --bins_obs 8 --N_iter 100 --no_zeros
```
The --no_zeros flag tells the code not to consider objects that do not satisfy any of the probability thresholds. If you don't explicitely state this, objects with p_cs < 0.8 and p_el < 0.8 will be considered as a third class.

The output of the run defined above should be something like:
```
L =  0.142363503568  +-  0.00203178319883
```
This shows the mean and standard deviation of the bias over 100 bootstrapping samples.

I am aware more documentation is needed and I'm currently working on it. I would be happy to assist you on how to run the code, so please contact me if you have some data that you would like to know how biased it is.
