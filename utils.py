import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

save_files = False

def default_val(key, col, default):
    if(key in col):
        return col[key]
    else:
        return default

#Takes a variable list of tuples in the form of ("Name", x_vals, y_vals) and plots them together on a graph
def plot(args, **kwargs):
    fig = plt.figure()

    fig_plot = fig.add_subplot(111)

    plots = []
    names = []

    colors = ["b", "g", "r", "c", "m", "y", "k"]
    for i, arg in enumerate(args):
        
        p, = fig_plot.plot(arg[1], arg[2], colors[i % 7], label=arg[0], marker=default_val("marker", kwargs, 'x'), linestyle=default_val("linestyle", kwargs, '-'))

        plots.append(p)
        names.append(arg[0])

    fig_plot.legend(plots, names, loc=2)

    fig_plot.set_xlabel(default_val("xlabel", kwargs, "x"))
    fig_plot.set_ylabel(default_val("ylabel", kwargs, "y"))
    fig_plot.set_title(default_val("title", kwargs, " vs. ".join(names)))
    
    fig.show()
    if save_files and "filename" in kwargs:
        plt.savefig(kwargs["filename"])
    

    raw_input("press enter to continue")
