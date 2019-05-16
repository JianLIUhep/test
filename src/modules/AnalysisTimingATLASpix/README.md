# AnalysisTimingATLASpix
**Maintainer**: Jens Kroeger (jens.kroeger@cern.ch)  
**Module Type**: *DUT*  
**Detector Type**: *ATLASpix*  
**Status**: work in progress

### Description
This module contains everything that's necessary for an in-depth timing analysis of the ATLASpix, including row and timewalk corrections.

Before being able to apply the row and timewalk correction, correction files need to be provided.
These can be generated by setting `calc_corrections=true`.
The calculation of the timewalk correction is based on a row-corrected histogram such that this procedure needs to be performed in 2 steps.

1. Calculated row correction.
2. Use row correction file and calculated timewalk correction.
After this both corrections can be applied on top of each other.

### Parameters
* `timing_cut`: Timing cut for associating a track with an ATLASpix cluster. Defaults to `1us`.
* `chi2ndof_cut`: Acceptance criterion for telescope tracks, defaults to a value of `3`.
* `time_cut_frame_edge`: Parameter to discard telescope tracks at the frame edges (start and end of the current frame). Defaults to `20ns`.
* `cluster_charge_cut`: Parameter to discard clusters with a charge larger than the cut. Defaults to `100000e` (inifitely large).
* `cluster_size_cut`: Parameter to discard clusters with a size too large, only for debugging purposes, default is 100 (inifitely large).
* `high_tot_cut`: Cut dividing 'low' and 'high' ToT events (based on seed pixel ToT). Defaults to `40`.
* `high_charge_cut`: Cut dividing 'low' and 'high' charge events (based on cluster charge). Defaults to `40e`.
* `left_tail_cut`: Cut to divide into left tail and main peak of time correlation histogram. Only used to investigate characteristics of left tail. Defaults to `-10ns`.
* `calc_corrections`: If `true`, TGraphErrors for row and timewalk corrections are produced.
* `correction_file_row`, `correction_file_timewalk`. Defaults to `false`.
* `correction_file_row`: Path to file which contains TGraphErrors for row correction. If this parameter is set, also `correction_graph_row` needs to be set. No default.
* `correction_file_timewalk`: Path to file which contains TGraphErrors for timewalk correction. If this parameter is set, also `correction_graph_timewalk` needs to be set. No default.
* `correction_graph_row`: Name of the TGraphErrors including its path in the root file used for row correction. E.g. "AnalysisTimingATLASpix/apx_0/gTimeCorrelationVsRow". No default.
* `correction_graph_timewalk`: Name of the TGraphErrors including its path in the root file used for row correction. E.g. "AnalysisTimingATLASpix/apx_0/gTimeCorrelationVsTot". No default.

### Plots produced
* 1D histograms
  * Track time correlation (all clusters)
  * Track time correlation (track-associated clusters)
  * Track time correlation (after row correction)
  * Track time correlation (after row and timewalk correction)
  * Track time correlation (after row and timewalk correction) for clusterToT < 25 lsb
  * Track time correlation (after row and timewalk correction) for clusterToT < 40 lsb
  * Track time correlation (after row and timewalk correction) for clusterToT > 40 lsb
  * Track time correlation example slice of 2D plot to investigate quality of Gaussian fit
  * Cluster time minus pixel time (for all pixels in cluster) to confirm that cluster time = time of earliest pixel
  * Pixel ToT for left tail events in time correlation (track timestamp - cluster timestamp < left_tail_cut)
  * Pixel ToT for main peak events in time correlation (track timestamp - cluster timestamp > left_tail_cut)
  * Pixel timestamp for left tail events in time correlation (track timestamp - cluster timestamp < left_tail_cut)
  * Pixel timestamp main peak events in time correlation (track timestamp - cluster timestamp > left_tail_cut)
  * Cluster size for left tail events in time correlation (track timestamp - cluster timestamp < left_tail_cut)
  * Cluster size for main peak events in time correlation (track timestamp - cluster timestamp > left_tail_cut)

* 2D histograms
  * Track time correlation vs. column (only control plot)
  * Track time correlation vs. cluster row (for row correction)
  * Track time correlation vs. cluster row (only single pixel clusters)
  * Track time correlation vs. cluster row (only multi-pixel clusters)
  * Track time correlation vs. cluster row (after row correction)
  * Track time correlation vs. seed pixel ToT (for timewalk correction)
  * Track time correlation vs. seed pixel ToT (only single pixel clusters)
  * Track time correlation vs. seed pixel ToT (only multi-pixel clusters)
  * Track time correlation vs. cluster row (after row correction)
  * Track time correlation vs. cluster row (after row correction, only single pixel clusters)
  * Track time correlation vs. cluster row (after row correction, only multi-pixel clusters)
  * Track time correlation vs. cluster row (after row and timewalk correction)
  * Track time correlation vs. seed pixel (after row and timewalk correction)
  * Cluster size vs. cluster ToT (only associated clusters)
  * Hit map of all pixels from associated clusters
  * Hit map of all pixels from associated clusters with high ToT
  * In-pixel distribution of tracks
  * In-pixel distribution of tracks (for clusters with high ToT)
  * Map of associated clusters
  * Pixel ToT vs. time for all clusters
  * Pixel ToT vs. time for high ToT clusters
  * Cluster map for left tail events in time correlation (track timestamp - cluster timestamp < left_tail_cut)
  * Cluster map for main peak events in time correlation (track timestamp - cluster timestamp > left_tail_cut)

* TGraphErrors
  * Peak of time correlation vs. row
  * Peak of time correlation vs. ToT (after row correction)
  * Peak of time correlation vs. row (after row correction, only single pixel clusters)
  * Peak of time correlation vs. row (after row correction, only multi-pixel clusters)

### Usage
```toml
[AnalysisTiming]
calc_corrections = false
correction_file_row = "correction_files/row_correction_file.root"
correction_graph_row = "AnalysisTimingATLASpix/apx0/gRTimeCorrelationVsRow"
correction_file_timewalk = "correction_files/timewalk_correction_file.root"
correction_graph_timewalk = "AnalysisTimingATLASpix/apx0/gTimCorrelationVsTot"
```