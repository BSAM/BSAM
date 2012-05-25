import matplotlib.pyplot as plt
import numpy
import bsam_io
import bsam_plot

# Parse datafile
data = bsam_io.parse_adaptive_data_file("out/m00001.dat")
root_patch = data['levels'][0][0]

# Create plots
fig = plt.figure()

ax  = fig.add_subplot(111)
levels = bsam_plot.create_levels(-0.1,0.7,50)
cs1 = bsam_plot.plot_contour(ax, root_patch, 3, levels)
fig.colorbar(cs1)
bsam_plot.plot_bounding_boxes(ax, data)

plt.show()
