import matplotlib.pyplot as plt
import numpy
import bsam_io
import bsam_plot

# Parse datafile
data = bsam_io.parse_adaptive_data_file("out/m00001.dat")
root_patch = data['levels'][0][0]

# Create plots
fig = plt.figure(figsize=(18,14))

ax  = fig.add_subplot(221)
levels = bsam_plot.create_levels(-0.1,0.7,50)
bsam_plot.plot_contourf(ax, root_patch, 0, levels)

ax  = fig.add_subplot(222)
cs1 = bsam_plot.plot_contourf(ax, root_patch, 3, levels)
fig.colorbar(cs1)
#bsam_plot.plot_bounding_boxes(ax, data)

ax  = fig.add_subplot(223)
bsam_plot.plot_xprofile(ax, root_patch,
                        [(0, 'k'),
                         (3, 'k--')],
                        [(-0.05,0.3),
                         (-1.0,1.0)])

ax  = fig.add_subplot(224)
bsam_plot.plot_xprofile(ax, root_patch,
                        [(0, 'k'),
                         (3, 'k--')])

plt.show()
