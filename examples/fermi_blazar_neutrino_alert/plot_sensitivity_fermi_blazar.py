import argparse
import matplotlib.pyplot as plt
import os

from alertstack.stats import TSHandler
from examples.fermi_blazar_neutrino_alert import blazar_analysis

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='Plot TS distribution and sensitivity'
    )
    parser.add_argument(
        '--input',
        type=str,
        default=(
            "/data/user/gsommani/alertstack-icecube/examples/"
            "fermi_blazar_neutrino_alert/cache/"
            "february_update_2026_02_16-09_37_22.pkl"
        ),
        help = 'Results to use')
    args = parser.parse_args()
    '''
    input: input file with the results to use.
    '''

    cwd = os.path.dirname(os.path.realpath(__file__))
    figures_folder = os.path.join(cwd,"figures/")
    
    res = blazar_analysis.load_results(filename=args.input)
    ts_handler = TSHandler(res, blazar_analysis)
    gd = ts_handler.find_thresholds_gamma()
    key = list(ts_handler.sens_threshold.keys())[0]
    val = res[0][key]
    ts_handler.plot_ts(val, key, gd=gd, bins=30) 
    plt.title(f"4LAC-DR3 (12y integrated flux) + IceCat-2 -> {len(
        val
    ):.1e} Scrambles")
    plt.savefig(
        figures_folder + "4LACDR3_integrated_scrambles_update_signalness",
        bbox_inches="tight",
        dpi=200
    )
    plt.savefig(
        figures_folder + "4LACDR3_integrated_scrambles_update_signalness.pdf",
        bbox_inches="tight",
        dpi=200
    )
    plt.close()
    ts_handler.extract_sens_dp(extent=0.15)
    ts_handler.plot_sens_dp(data_derived=False)
    
    plt.legend(loc=(-0.6, 0))
    
    plt.text(
        ts_handler.x1 + 0.001, 0.95,
        f"{ts_handler.x1*100:.1f}%",
        color="tab:blue",
        rotation=90
    )
    plt.text(
        ts_handler.x2 + 0.001, -0.03,
        f"{ts_handler.x2*100:.1f}%",
        color="tab:orange",
        rotation=90
    )
    plt.text(
        ts_handler.x3 - 0.005, -0.03,
        f"{ts_handler.x3*100:.1f}%",
        color="tab:green",
        rotation=90
    )
    
    plt.text(-0.005, 0.85, "90%", color="black")
    plt.text(-0.005, 0.45, "50%", color="black")
    
    plt.savefig(
        figures_folder + "4LACDR3_integrated_sensitivity_5sigma_update_signalness",
        bbox_inches="tight",
        dpi=200
    )
    plt.savefig(
        figures_folder + "4LACDR3_integrated_sensitivity_5sigma_update_signalness.pdf",
        bbox_inches="tight",
        dpi=200
    )
    plt.close()
