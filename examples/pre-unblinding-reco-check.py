import matplotlib.pyplot as plt
import numpy as np
import os

from alertstack.stats import TSHandler
from examples.fermi_blazar_neutrino_alert import (
    blazar_analysis as blazar_analysis_average
)
from examples.fermi_LC_blazar_neutrino_alert import (
    blazar_analysis as blazar_analysis_monthly
)
from examples.radio_agn_neutrino_alert import agn_analysis
from matplotlib.ticker import AutoMinorLocator


result_fermi_average_old_hed = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/fermi_blazar_neutrino_alert/cache/"
    "2026_03_31-13_05_28.pkl"
)
result_fermi_monthly_old_hed = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/fermi_LC_blazar_neutrino_alert/cache/"
    "2026_03_31-14_03_47.pkl"
)
result_radio_old_hed = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/radio_agn_neutrino_alert/cache/"
    "2026_03_31-11_54_05.pkl"
)
result_fermi_average_old_led = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/fermi_blazar_neutrino_alert/cache/"
    "2026_03_30-14_52_32.pkl"
)
result_fermi_monthly_old_led = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/fermi_LC_blazar_neutrino_alert/cache/"
    "2026_03_30-18_44_43.pkl"
)
result_radio_old_led = (
    "/data/user/gsommani/wg-nu-sources/2022_Fermi_Blazars_Alerts_Stacking/"
    "examples/radio_agn_neutrino_alert/cache/"
    "2026_03_31-04_11_11.pkl"
)


result_fermi_average_new_hed = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/fermi_blazar_neutrino_alert/cache/"
    "HED_2026_03_31-08_22_00.pkl"
)
result_fermi_monthly_new_hed = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/fermi_LC_blazar_neutrino_alert/cache/"
    "HED_2026_03_31-08_30_50.pkl"
)
result_radio_new_hed = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/radio_agn_neutrino_alert/cache/"
    "HED_2026_03_31-08_23_32.pkl"
)
result_fermi_average_new_led = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/fermi_blazar_neutrino_alert/cache/"
    "LED_2026_03_31-08_13_02.pkl"
)
result_fermi_monthly_new_led = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/fermi_LC_blazar_neutrino_alert/cache/"
    "LED_2026_03_31-08_23_27.pkl"
)
result_radio_new_led = (
    "/data/user/gsommani/alertstack-icecube/"
    "examples/radio_agn_neutrino_alert/cache/"
    "LED_2026_03_31-08_13_41.pkl"
)

def extract_sensitivity(
    analysis,
    result_path,
    max_run=134818,
    evttype="ALL",
    extent=0.4,
    only_sens=True
):
    """Given the results from an analysis, extract the sensitivity
    in terms of percentage of diffuse astrophysical flux.

    Parameters
    ----------
    analysis: `alertstack.analysis.Analyse`
        The analysis to which the results refer
    result_path: `str`
        Path to the results to use
    max_run: `int`
        Last run to consider
    evttype: `str`
        Select all neutrinos ('ALL'), only LED neutrinos ('LED),
        or only HED neutrinos ('HED')
    extent: `float`
        Extent to use for interpolations (in terms of astrophysical
        neutrino flux)
    only_sens: `bool`
        If True, only the sensitivity is extrapolated and not the
        discovery potentials.
    """
    result = analysis.load_results(filename=result_path)
    ts_handler = TSHandler(
        result, analysis, max_run=max_run, evttype=evttype
    )
    ts_handler.find_thresholds_gamma()
    ts_handler.extract_sens_dp(extent=extent, only_sens=only_sens)
    sensitivity = ts_handler.x1
    return sensitivity

if __name__ == "__main__":

    cwd = os.path.dirname(os.path.realpath(__file__))
    figures_folder = os.path.join(
        cwd,"fermi_blazar_neutrino_alert/figures/"
    )

    
    sens_fermi_average_old_hed = extract_sensitivity(
        blazar_analysis_average,
        result_fermi_average_old_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Average Old HED: {sens_fermi_average_old_hed}\n\n")
    sens_fermi_monthly_old_hed = extract_sensitivity(
        blazar_analysis_monthly,
        result_fermi_monthly_old_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Monthly Old HED: {sens_fermi_monthly_old_hed}\n\n")
    sens_radio_old_hed = extract_sensitivity(
        agn_analysis,
        result_radio_old_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tRadio Old HED: {sens_radio_old_hed}\n\n")
    sens_fermi_average_old_led = extract_sensitivity(
        blazar_analysis_average, result_fermi_average_old_led, evttype="LED"
    )
    print(f"\n\n\tFermi Average Old LED: {sens_fermi_average_old_led}\n\n")
    sens_fermi_monthly_old_led = extract_sensitivity(
        blazar_analysis_monthly, result_fermi_monthly_old_led, evttype="LED"
    )
    print(f"\n\n\tFermi Monthly Old LED: {sens_fermi_monthly_old_led}\n\n")
    sens_radio_old_led = extract_sensitivity(
        agn_analysis, result_radio_old_led, evttype="LED"
    )
    print(f"\n\n\tRadio Old LED: {sens_radio_old_led}\n\n")

    
    sens_fermi_average_new_hed = extract_sensitivity(
        blazar_analysis_average,
        result_fermi_average_new_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Average New HED: {sens_fermi_average_new_hed}\n\n")
    sens_fermi_monthly_new_hed = extract_sensitivity(
        blazar_analysis_monthly,
        result_fermi_monthly_new_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Monthly New HED: {sens_fermi_monthly_new_hed}\n\n")
    sens_radio_new_hed = extract_sensitivity(
        agn_analysis,
        result_radio_new_hed,
        evttype="HED",
        extent=0.15,
    )
    print(f"\n\n\tRadio New HED: {sens_radio_new_hed}\n\n")
    sens_fermi_average_new_led = extract_sensitivity(
        blazar_analysis_average,
        result_fermi_average_new_led,
        evttype="LED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Average New LED: {sens_fermi_average_new_led}\n\n")
    sens_fermi_monthly_new_led = extract_sensitivity(
        blazar_analysis_monthly,
        result_fermi_monthly_new_led,
        evttype="LED",
        extent=0.15,
    )
    print(f"\n\n\tFermi Monthly New LED: {sens_fermi_monthly_new_led}\n\n")
    sens_radio_new_led = extract_sensitivity(
        agn_analysis,
        result_radio_new_led,
        evttype="LED",
        extent=0.15,
    )
    print(f"\n\n\tRadio New LED: {sens_radio_new_led}\n\n")

    impr_fermi_average_led = (
        sens_fermi_average_old_led / sens_fermi_average_new_led
    )
    impr_fermi_average_hed = (
        sens_fermi_average_old_hed / sens_fermi_average_new_hed
    )

    impr_fermi_monthly_led = (
        sens_fermi_monthly_old_led / sens_fermi_monthly_new_led
    )
    impr_fermi_monthly_hed = (
        sens_fermi_monthly_old_hed / sens_fermi_monthly_new_hed
    )

    impr_radio_led = (
        sens_radio_old_led / sens_radio_new_led
    )
    impr_radio_hed = (
        sens_radio_old_hed / sens_radio_new_hed
    )
    
    print("\n\nFermi Average")
    print(f"\tImprovement LED: {impr_fermi_average_led}")
    print(f"\tImprovement HED: {impr_fermi_average_hed}")
    
    print("\n\nFermi Monthly")
    print(f"\tImprovement LED: {impr_fermi_monthly_led}")
    print(f"\tImprovement HED: {impr_fermi_monthly_hed}")
    
    print("\n\nRadio")
    print(f"\tImprovement LED: {impr_radio_led}")
    print(f"\tImprovement HED: {impr_radio_hed}")
    print("\n\n")

    
    xs = [1,2,3]
    plt.xlim(0.5,3.5)
    plt.ylim(0.5,2.5)
    #plt.grid(axis="x")
    plt.grid(linestyle="dashed", axis="y", which="major", linewidth=0.5)
    plt.grid(linestyle="dotted", axis="y", which="minor", linewidth=0.5)
    plt.scatter(xs, [
        impr_radio_led,
        impr_fermi_average_led,
        impr_fermi_monthly_led,
    ], label="SplineMPE with likelihood scan")
    plt.scatter(xs, [
        impr_radio_hed,
        impr_fermi_average_hed,
        impr_fermi_monthly_hed,
        ], label="Millipede Wilks"
    )

    plt.axhline(1., color="black", label="No improvement")

    axs = plt.gca()
    axs.yaxis.set_minor_locator(AutoMinorLocator())
    axs.set_xticks(xs, labels=["RFC", "Fermi average", "Fermi monthly"])
    axs.set_axisbelow(True)
    plt.ylabel("Old / New sensitivity")
    plt.legend()
    
    plt.savefig(
        os.path.join(figures_folder,"led_hed_sens_comparison"),
        dpi=150,
        bbox_inches='tight'
    )
    plt.close()

    
    xs = np.array([1,2,3])
    xs_led = xs - 0.0
    xs_hed = xs + 0.0
    plt.xlim(0.5,3.5)
    plt.ylim(5,30)
    #plt.grid(axis="x")
    plt.grid(linestyle="dashed", axis="y", which="major", linewidth=0.5)
    plt.grid(linestyle="dotted", axis="y", which="minor", linewidth=0.5)
    plt.scatter(xs_led, np.array([
        sens_radio_old_led,
        sens_fermi_average_old_led,
        sens_fermi_monthly_old_led,
    ])*100, label="Cristina's analysis, LED alerts"
    )
    plt.scatter(xs_hed, np.array([
        sens_radio_old_hed,
        sens_fermi_average_old_hed,
        sens_fermi_monthly_old_hed,
        ])*100, label="Cristina's analysis, HED alerts"
    )
    plt.scatter(xs_led, np.array([
        sens_radio_new_led,
        sens_fermi_average_new_led,
        sens_fermi_monthly_new_led,
    ])*100, label="This analysis, LED alerts", color="tab:blue", marker="*"
    )
    plt.scatter(xs_hed, np.array([
        sens_radio_new_hed,
        sens_fermi_average_new_hed,
        sens_fermi_monthly_new_hed,
        ])*100, label="This analysis, HED alerts", color="tab:orange", marker="*"
    )

    axs = plt.gca()
    plt.yscale('log')
    #axs.yaxis.set_minor_locator(AutoMinorLocator())
    axs.set_xticks(xs, labels=["RFC", "Fermi average", "Fermi monthly"])
    plt.ylabel("Percentage of astrophysical neutrino flux [%]")
    plt.legend()
    
    plt.savefig(
        os.path.join(figures_folder,"led_hed_sens_comparison_abs"),
        dpi=150,
        bbox_inches='tight'
    )
    plt.close()
