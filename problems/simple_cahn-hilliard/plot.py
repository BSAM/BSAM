import matplotlib.pyplot as plt
import numpy
import bsam_io
import bsam_plot

# Create plots
fig = plt.figure()

ax  = fig.add_subplot(211)
data = bsam_io.parse_adaptive_data_file("out/m00000.dat")
root_patch = data['levels'][0][0]
bsam_plot.plot_grid(ax, data, 5)
bsam_plot.plot_interface(ax, root_patch)
bsam_plot.plot_contourf(ax, root_patch, 0)

ax  = fig.add_subplot(212)
data = bsam_io.parse_adaptive_data_file("out/m00010.dat")
root_patch = data['levels'][0][0]
#bsam_plot.plot_grid(ax, data, 5)
#bsam_plot.plot_interface(ax, root_patch)
bsam_plot.plot_contourf(ax, root_patch, 0)

plt.show()
