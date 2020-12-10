# Guenzel-McCollum-et-al
Analysis scripts for Günzel & McCollum et al.; Social modulation of individual preferences in cockroaches; iScience
## PreferenceInversionInCockroaches_Behaviour.m.
This script represents the pipeline used to analyse behavioural data. It first pools data across different trials, then analyses them in terms of, for example, how much time animals spent in any of the two shelters, last the script creates figures and calculates statistics.
## PreferenceInversionInCockroaches_Physiology.m
This script represents the pipeline used to analyse calcium-imaging data. It uses baseline-subtracted Delta R340/380 2D fluorescence time courses and applies a k-medoids cluster analysis on active pixels (fluorescence > 0.01 for at least one of the five odours).
