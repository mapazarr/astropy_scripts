def plot_something(x, y, ax=None, style_kwargs=None):
    """Plot something.
    
    Parameters
    ----------
    x, y : tbd
        tbd
    ax : `~matplotlib.axes.Axes` or None
        matplotlib axes to plot on
    style_kwargs : dict or None
        Style options, passed to matplotlib
    
    Returns
    -------
    ax : `~matplotlib.axes.Axes` or None
        matplotlib axes to plot on
    """
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()
    if style_kwargs is None:
        style_kwargs = dict()

    # Plotting happens here
    ax.plot(x, y, **style_kwargs)

    return ax
