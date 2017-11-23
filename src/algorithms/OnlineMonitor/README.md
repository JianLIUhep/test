## OnlineMonitor
**Maintainer**: Daniel Hynds (<daniel.hynds@cern.ch>)   
**Status**: Functional - some minor fixes needed

#### Description
This algorithm opens a GUI to monitor the progress of Corryvreckan reconstruction. Each canvas contains real time plots of the reconstruction, updated every 500 events. Plots are used from the `BasicTracking` and `TestAlgorithm` algorithms.

Users should be able to exit `OnlineMonitor` and leave the reconstruction still running, but currently this causes a crash.

#### Parameters
* `reference`: Name of the reference plane.

#### Plots produced
Overview canvas:
* Cluster ToT of reference plane
* 2D hitmap of reference plane
* Residual in X of reference plane

Tracking canvas:
* Track chi^2
* Track angle in X

For each detector the following plots are produced:
* Hitmap canvas: 2D hitmap
* Residuals canvas: residual in X histogram
* Event times canvas: event times histogram
* Correlations X canvas: correlation in X plot
* Correlations Y canvas: correlation in Y plot
* Charge distributions canvas: cluster ToT - broken at the moment

#### Usage
```toml
[OnlineMonitor]
reference = "W0013_E03"
```
Parameters to be used in multiple algorithms can also be defined globally at the top of the configuration file. This is highly encouraged for parameters such as `DUT` and `reference`.