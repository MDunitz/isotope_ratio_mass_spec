import bokeh
from bokeh.models import Panel, HoverTool, LabelSet, ColumnDataSource, Tabs

from bokeh.plotting import figure, output_file, save
import time


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
        tools=[HoverTool(tooltips=[("Compound", "@compound")])],
    )
    p.circle(source=df, x="relative_time", y="original_peak_order", color="color")
    return p


def plot_relative_time_V_peak_order_by_tab(df):

    titles = list(df["Analysis"].unique())

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
        width=1000,
        height=600,
        x_axis_label="Sample",
        y_axis_label="Peak Count",
        toolbar_location="above",
    )
    p.circle(source=df, x="bio_replicate_id", y="original_peak_count", color="color")
    return p


def plot_relative_time_V_amplitude_with_labels(df, prelabeled):
    source = ColumnDataSource(
        df[["relative_time", "amplitude_2", "TENTATIVE_COMPOUND", "color"]]
    )
    control_source = ColumnDataSource(
        prelabeled[["relative_time", "amplitude_2", "compound", "color"]]
    )

    p = bokeh.plotting.figure(
        width=1200,
        height=600,
        x_axis_label="Retention Time",
        y_axis_label="Amplitude",
        toolbar_location="above",
        # tools=[
        #     HoverTool(
        #         names=["tent_comp"], tooltips=[("Compound", "@TENTATIVE_COMPOUND")]
        #     ),
        #     HoverTool(names=["prelabel"], tooltips=[("Prelabeled", "@compound")]),
        # ],
    )
    # points with names created in source

    p.vbar(
        source=source,
        x="relative_time",
        top="amplitude_2",
        width=0.5,
        bottom=0,
        color="blue",
    )
    p.circle(
        source=source,
        x="relative_time",
        y="amplitude_2",
        color="blue",
        name="tent_comp",
    )

    labels = LabelSet(
        x="relative_time",
        y="amplitude_2",
        text="TENTATIVE_COMPOUND",
        x_offset=3,
        y_offset=3,
        source=source,
        text_font_size={"value": "16px"},
    )

    # Adding that label to our figure
    p.add_layout(labels)

    # include prelabeled set
    p.square(
        source=control_source,
        x="relative_time",
        y="amplitude_2",
        color="red",
        name="prelabel",
    )
    control_labels = LabelSet(
        x="relative_time",
        y="amplitude_2",
        text="compound",
        x_offset=3,
        y_offset=3,
        source=control_source,
        text_font_size={"value": "16px"},
        text_color="red",
    )
    p.add_layout(control_labels)
    return p


def plot_relative_time_V_peak_amplitude_by_tab_with_labels(
    df, save=False, filename=str(time.time())
):
    print(filename)
    titles = list(df["Analysis"].unique())
    prelabeled = df[df["Analysis"] == 4167]
    fig_list = []
    for x in titles:
        subset = df[df["Analysis"] == x]
        fig_list.append(plot_relative_time_V_amplitude_with_labels(subset, prelabeled))

    # Generating tab image from list
    tabs = [
        Panel(child=fig_list[l], title=str(titles[l])) for l in range(len(fig_list))
    ]
    if save:
        save_html_to_figures(tabs, filename)
    return tabs


def save_html_to_figures(p, title):
    output_file(f"figures/{title}.html")
    save(Tabs(tabs=p))
