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

cristina_tot_nu = 275
cristina_average_signalness = 0.451
cristina_sens_nu_fermi_average = 6.3
cristina_sens_nu_fermi_monthly = 5.5
cristina_sens_nu_rfc = 7.4
result_fermi_average = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "fermi_blazar_neutrino_alert/cache/"
    "february_update_2026_02_27-11_49_46.pkl"
)
result_fermi_monthly = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "fermi_LC_blazar_neutrino_alert/cache/"
    "february_monthly_2026_02_27-13_02_24.pkl"
)
result_agn = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "radio_agn_neutrino_alert/cache/"
    "february_update_2026_02_27-12_01_55.pkl"
)
result_fermi_average_reduced = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "fermi_blazar_neutrino_alert/cache/"
    "reducedalerts_2026_03_20-11_02_29.pkl"
)
result_fermi_monthly_reduced = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "fermi_LC_blazar_neutrino_alert/cache/"
    "reducedalerts_2026_03_20-10_50_35.pkl"
)
result_agn_reduced = (
    "/data/user/gsommani/alertstack-icecube/examples/"
    "radio_agn_neutrino_alert/cache/"
    "reducedalerts_2026_03_20-11_04_48.pkl"
)

if __name__ == "__main__":

    cwd = os.path.dirname(os.path.realpath(__file__))
    figures_folder = os.path.join(
        cwd,"fermi_blazar_neutrino_alert/figures/"
    )

    res_average = blazar_analysis_average.load_results(
        filename=result_fermi_average
    )
    res_monthly = blazar_analysis_monthly.load_results(
        filename=result_fermi_monthly
    )
    res_rfc = agn_analysis.load_results(
        filename=result_agn
    )
    res_average_reduced = blazar_analysis_average.load_results(
        filename=result_fermi_average_reduced
    )
    res_monthly_reduced = blazar_analysis_monthly.load_results(
        filename=result_fermi_monthly_reduced
    )
    res_rfc_reduced = agn_analysis.load_results(
        filename=result_agn_reduced
    )
    ts_handler_average = TSHandler(res_average, blazar_analysis_average)
    ts_handler_monthly = TSHandler(res_monthly, blazar_analysis_monthly)
    ts_handler_rfc = TSHandler(res_rfc, agn_analysis)
    ts_handler_average_reduced = TSHandler(
        res_average_reduced, blazar_analysis_average
    )
    ts_handler_monthly_reduced = TSHandler(
        res_monthly_reduced, blazar_analysis_monthly
    )
    ts_handler_rfc_reduced = TSHandler(res_rfc_reduced, agn_analysis)
    ts_handler_average.find_thresholds_gamma()
    ts_handler_monthly.find_thresholds_gamma()
    ts_handler_rfc.find_thresholds_gamma()
    ts_handler_average_reduced.find_thresholds_gamma()
    ts_handler_monthly_reduced.find_thresholds_gamma()
    ts_handler_rfc_reduced.find_thresholds_gamma()
    ts_handler_average.extract_sens_dp(extent=0.15)
    ts_handler_monthly.extract_sens_dp(extent=0.15)
    ts_handler_rfc.extract_sens_dp(extent=0.15)
    ts_handler_average_reduced.extract_sens_dp(extent=0.175)
    ts_handler_monthly_reduced.extract_sens_dp(extent=0.175)
    ts_handler_rfc_reduced.extract_sens_dp(extent=0.2)
    sens_average = ts_handler_average.x1
    sens_monthly = ts_handler_monthly.x1
    sens_rfc = ts_handler_rfc.x1
    sens_average_reduced = ts_handler_average_reduced.x1
    sens_monthly_reduced = ts_handler_monthly_reduced.x1
    sens_rfc_reduced = ts_handler_rfc_reduced.x1

    cristina_tot_signal_nu = cristina_tot_nu * cristina_average_signalness
    cristina_sens_fermi_average = (
        cristina_sens_nu_fermi_average / cristina_tot_signal_nu
    )
    cristina_sens_fermi_monthly = (
        cristina_sens_nu_fermi_monthly / cristina_tot_signal_nu
    )
    cristina_sens_rfc = (
        cristina_sens_nu_rfc / cristina_tot_signal_nu
    )

    xs = [1,2,3]
    plt.xlim(0.5,3.5)
    plt.ylim(0,7)
    plt.grid(axis="x")
    plt.grid(linestyle="dotted", axis="y")
    plt.scatter(xs, np.array([
        cristina_sens_rfc,
        cristina_sens_fermi_average,
        cristina_sens_fermi_monthly,
    ])*100, label="Cristina's sensitivities")
    plt.scatter(
        xs,
        np.array(
            [sens_rfc_reduced, sens_average_reduced, sens_monthly_reduced]
        )*100,
        label="New recos"
    )
    plt.scatter(
        xs,
        np.array([sens_rfc, sens_average, sens_monthly])*100,
        label="New recos + new alerts"
    )

    axs = plt.gca()
    axs.set_xticks(xs, labels=["RFC", "Fermi average", "Fermi monthly"])
    plt.ylabel("Percentage of astrophysical neutrino flux [%]")
    plt.legend()
    
    plt.savefig(
        os.path.join(figures_folder,"cristina_sens_comparison"),
        dpi=150,
        bbox_inches='tight'
    )
    plt.close()
    