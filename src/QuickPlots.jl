#=

    Functions for quick plotting using CairoMakie with pre-defined settings.

=#

module QuickPlots
using CairoMakie
using ..PlotSettings

# plot settings 
CairoMakie.activate!(;px_per_unit=4)
PlotSettings.set_global_themes!()
text_xloc = 0.05
text_yloc = 0.95
xalign = :left
yalign = :top;
fig_width_1, fig_height_1, fig_width_1_2, fig_height_1_2, fig_width_2_2, fig_height_2_2, fontsize, xlabelsize, ylabelsize, xticklabelsize, yticklabelsize, col_gap, xgridvisible, ygridvisible = PlotSettings.load_settings();

default_colors = [:tomato, :royalblue, :aquamarine4, :black, :rebeccapurple, :lightpink2]
default_labels = ["", "", "", "", "", "", ""]
default_linestyles = [:solid, :solid, :solid, :solid, :solid, :solid]
default_linewidths = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
default_alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

alt_linestyles = [:solid, :dash, :dashdot, :dashdotdot, :dot, :dash]
alt_linewidths = [2.0, 1.9, 1.8, 1.7, 1.6, 1.5]
alt_alphas = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]

function plot11(x, y;
    colors = default_colors,
    labels = default_labels,
    linestyles = default_linestyles,
    linewidths = default_linewidths,
    alphas = default_alphas,
    xlabel = "",
    ylabel = "",
    lim_x_min = nothing,
    lim_x_max = nothing,
    lim_y_min = nothing,
    lim_y_max = nothing,
    xscale = identity,
    yscale = identity,
    text_xloc = nothing, # loacation of text overlayed on plot
    text_yloc = nothing, # loacation of text overlayed on plot
    text = nothing, # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :rt, # position of legend
    labelsize = fontsize, # size of legend and axis labels
    framevisible = true, # whether to draw box around legend
    hlines = nothing, # y-values to draw horizontal lines at in the form [y1, y2, ...]
    vlines = nothing, # x-values to draw vertical lines at in the form [x1, x2, ...]
    save_plot = false, # whether to save the plot to file
    fname = "", # file name to save the plot to
    scatter_lines = false # whether to plot scatter lines instead of normal lines
    )

    with_theme(theme_latexfonts()) do
    f = Figure(size = (fig_width_1, fig_height_1));
    ax = Axis(f[1,1], 
        limits = (lim_x_min, lim_x_max, lim_y_min, lim_y_max),
        xscale = xscale,
        yscale = yscale,
        ylabel=ylabel,
        xlabel=xlabel,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,
        # xticks = xticks,
        # yticks = yticks,
    )


    for i in eachindex(x)
        if scatter_lines
            scatterlines!(ax, x[i], y[i], color=colors[i], label=labels[i])
        else
            lines!(ax, x[i], y[i], color=colors[i], label=labels[i], linestyle=linestyles[i], linewidth=linewidths[i], alpha=alphas[i])
        end
    end

    if typeof(text) <: Nothing
        nothing
    else
        text!(ax, text_xloc, text_yloc; text = text, space=:relative, align = (xalign, yalign), fontsize = fontsize)
    end

    if vlines !== nothing
        for vv in vlines
            vlines!(ax, [vv], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if hlines !== nothing
        for hh in hlines
            hlines!(ax, [hh], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if legend
        axislegend(ax, position = position, labelsize = labelsize, framevisible = framevisible)
    end

    resize_to_layout!(f)
    if save_plot
        save(fname, f)
    end
    display(f)
    end
end

function plot12(x_left, y_left, x_right, y_right;
    colors_left = default_colors,
    labels_left = default_labels,
    linestyles_left = default_linestyles,
    linewidths_left = default_linewidths,
    alphas_left = default_alphas,
    colors_right = default_colors,
    labels_right = default_labels,
    linestyles_right = default_linestyles,
    linewidths_right = default_linewidths,
    alphas_right = default_alphas,
    xlabel_left = "",
    ylabel_left = "",
    xlabel_right = "",
    ylabel_right = "",
    lim_x_min_left = nothing,
    lim_x_max_left = nothing,
    lim_y_min_left = nothing,
    lim_y_max_left = nothing,
    lim_x_min_right = nothing,
    lim_x_max_right = nothing,
    lim_y_min_right = nothing,
    lim_y_max_right = nothing,
    xscale_left = identity,
    yscale_left = identity,
    xscale_right = identity,
    yscale_right = identity,
    text_left = nothing,
    text_right = nothing,
    text_xloc_left = nothing,
    text_yloc_left = nothing,
    text_xloc_right = nothing,
    text_yloc_right = nothing,
    xalign_left = :left,
    yalign_left = :top,
    xalign_right = :left,
    yalign_right = :top,
    legend_left = false,
    legend_right = false,
    position_left = :rt,
    position_right = :rt,
    labelsize_left = fontsize,
    labelsize_right = fontsize,
    framevisible_left = true,
    framevisible_right = true,
    hlines_left = nothing,
    vlines_left = nothing,
    hlines_right = nothing,
    vlines_right = nothing,
    save_plot = false,
    fname = "",
    scatter_lines_left = false,
    scatter_lines_right = false,
    )

    with_theme(theme_latexfonts()) do
    f = Figure(size = (fig_width_1_2, fig_height_1_2));

    ax_left = Axis(f[1,1],
        limits = (lim_x_min_left, lim_x_max_left, lim_y_min_left, lim_y_max_left),
        xscale = xscale_left,
        yscale = yscale_left,
        ylabel = ylabel_left,
        xlabel = xlabel_left,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,
    )

    ax_right = Axis(f[1,2],
        limits = (lim_x_min_right, lim_x_max_right, lim_y_min_right, lim_y_max_right),
        xscale = xscale_right,
        yscale = yscale_right,
        ylabel = ylabel_right,
        xlabel = xlabel_right,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,
    )

    colgap!(f.layout, col_gap)

    for i in eachindex(x_left)
        if scatter_lines_left
            scatterlines!(ax_left, x_left[i], y_left[i], color=colors_left[i], label=labels_left[i])
        else
            lines!(ax_left, x_left[i], y_left[i], color=colors_left[i], label=labels_left[i], linestyle=linestyles_left[i], linewidth=linewidths_left[i], alpha=alphas_left[i])
        end
    end

    for i in eachindex(x_right)
        if scatter_lines_right
            scatterlines!(ax_right, x_right[i], y_right[i], color=colors_right[i], label=labels_right[i])
        else
            lines!(ax_right, x_right[i], y_right[i], color=colors_right[i], label=labels_right[i], linestyle=linestyles_right[i], linewidth=linewidths_right[i], alpha=alphas_right[i])
        end
    end

    if !(typeof(text_left) <: Nothing)
        text!(ax_left, text_xloc_left, text_yloc_left; text = text_left, space=:relative, align = (xalign_left, yalign_left), fontsize = fontsize)
    end

    if !(typeof(text_right) <: Nothing)
        text!(ax_right, text_xloc_right, text_yloc_right; text = text_right, space=:relative, align = (xalign_right, yalign_right), fontsize = fontsize)
    end

    if vlines_left !== nothing
        for vv in vlines_left
            vlines!(ax_left, [vv], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if hlines_left !== nothing
        for hh in hlines_left
            hlines!(ax_left, [hh], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if vlines_right !== nothing
        for vv in vlines_right
            vlines!(ax_right, [vv], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if hlines_right !== nothing
        for hh in hlines_right
            hlines!(ax_right, [hh], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if legend_left
        axislegend(ax_left, position = position_left, labelsize = labelsize_left, framevisible = framevisible_left)
    end

    if legend_right
        axislegend(ax_right, position = position_right, labelsize = labelsize_right, framevisible = framevisible_right)
    end

    resize_to_layout!(f)
    if save_plot
        save(fname, f)
    end
    display(f)
    end
end


end
