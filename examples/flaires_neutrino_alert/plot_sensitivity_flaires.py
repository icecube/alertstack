import argparse
import matplotlib.pyplot as plt
import os

from alertstack.stats import TSHandler
from examples.flaires_neutrino_alert import flaires_analysis

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='Plot TS distribution and sensitivity'
    )
    parser.add_argument(
        '--input',
        type=str,
        default=(
            "/data/user/gsommani/alertstack_results/"
            "flaires_n2000000_f0.07_s10_update_signalness.pkl"
        ),
        help = 'Results to use')
    args = parser.parse_args()
    '''
    input: input file with the results to use.
    '''

    cwd = os.path.dirname(os.path.realpath(__file__))
    figures_folder = os.path.join(cwd,"figures/")
    print(figures_folder)
    
    res = flaires_analysis.load_results(filename=args.input)
    ts_handler = TSHandler(res, flaires_analysis)
    ts_handler.find_thresholds_from_data()
    key = list(ts_handler.sens_threshold.keys())[0]
    val = res[0][key]
    ts_handler.plot_ts(val, key, bins=20) 
    plt.title(f"Flaires + IceCat-2 -> {len(val):.1e} Scrambles")
    plt.savefig(
        figures_folder + "Flaires_scrambles_update_signalness",
        bbox_inches="tight",
        dpi=200
    )
    plt.savefig(
        figures_folder + "Flaires_scrambles_update_signalness.pdf",
        bbox_inches="tight",
        dpi=200
    )
    plt.close()
    ts_handler.extract_sens_dp(extent=0.07)
    ts_handler.plot_sens_dp()
    
    plt.legend()
    
    plt.text(
        ts_handler.x1 + 0.001, -0.03,
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
        ts_handler.x3 - 0.003, -0.03,
        f"{ts_handler.x3*100:.1f}%",
        color="tab:green",
        rotation=90
    )
    
    plt.text(-0.003, 0.85, "90%", color="black")
    plt.text(-0.003, 0.45, "50%", color="black")
    
    plt.savefig(
        figures_folder + "flaires_sensitivity_5sigma_update_signalness",
        bbox_inches="tight",
        dpi=200
    )
    plt.savefig(
        figures_folder + "flaires_sensitivity_5sigma_update_signalness.pdf",
        bbox_inches="tight",
        dpi=200
    )
    plt.close()
