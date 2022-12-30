import bokeh
from bokeh.transform import jitter
from bokeh.models import Panel, Tabs


def plot_peak_order_V_amplitude_2(df):
    p = bokeh.plotting.figure(
        width=800,
        height=300,
        x_axis_label="Peak Order",
        y_axis_label="Amplitude",
        toolbar_location="above",
    )
    p.circle(source=df, x="peak_order", y="amplitude_2", color="color")
    return p


def plot_time_relative_to_alanine_V_amplitude_2(df):
    p = bokeh.plotting.figure(
        width=800,
        height=300,
        x_axis_label="Time relative to alanine",
        y_axis_label="Amplitude",
        toolbar_location="above",
    )
    p.circle(source=df, x="TIME_RELATIVE_TO_ALANINE", y="amplitude_2", color="color")
    return p


def plot_relative_time_V_peak_order(df):
    p = bokeh.plotting.figure(
        width=800,
        height=300,
        x_axis_label="Relative Time",
        y_axis_label="Peak Order",
        toolbar_location="above",
    )
    p.circle(source=df, x="relative_time", y="original_peak_order", color="color")
    return p


def plot_relative_time_V_peak_order_by_tab(df):

    titles = list(df["bio_replicate_id"].unique())

    fig_list = []
    for x in titles:
        subset = df[df["bio_replicate_id"] == x]
        fig_list.append(plot_relative_time_V_peak_order(subset))

    # Generating tab image from list
    tabs = [
        Panel(child=fig_list[l], title=str(titles[l])) for l in range(len(fig_list))
    ]
    # bokeh.io.show(Tabs(tabs=tabs))
    return tabs


def plot_sample_V_peak_count(df):
    p = bokeh.plotting.figure(
        width=800,
        height=300,
        x_axis_label="Sample",
        y_axis_label="Peak Count",
        toolbar_location="above",
    )
    p.circle(source=df, x="bio_replicate_id", y="original_peak_count", color="color")
    return p
