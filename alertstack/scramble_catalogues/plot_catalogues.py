import matplotlib.pyplot as plt
import os

from agn_radio_catalogue import AstrogeoAGNCatalogue
from alertstack.fixed_catalogues.icecube_neutrino_alerts import (
    HealpixNeutrinoAlertCatalogue
)
from alertstack.plotting_tools import plot_catalogue
from blazar_catalogue import (
    Fermi4FGLBlazarCatalogue,
    AverageFluxWeightHypothesis,
)
from flaires_catalogue import (
    FlairesCatalogue,
    FluencebolHypothesis,
)
from sjoert_catalogue import (
    AccretionFlaresSjoertCatalogue,
    StrengthFluxWeightHypothesis,
)


if __name__ == "__main__":

    cwd = os.path.dirname(os.path.realpath(__file__))
    figures_folder = os.path.join(cwd, 'figures/')

    # Sjoert's catalogue
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        weighter=StrengthFluxWeightHypothesis(
            HealpixNeutrinoAlertCatalogue()
        ),
        tw=True
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_data_weight"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_data_weight.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_data"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_data.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        scramble=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_scramble_1"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        scramble=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_scramble_2"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        scramble=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_scramble_3"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        nside=4,
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_binning"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        AccretionFlaresSjoertCatalogue(),
        15,
        "Accretion Flares",
        "63 Accretion Flares",
        final_pdf=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"sjoert_accr_flares_pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()


    # Flaires catalogue
    catalogue_obj = FlairesCatalogue()
    plot_catalogue(
        catalogue_obj,
        5,
        "IR Flares",
        f"{len(catalogue_obj.parse_data())} sources",
        weighter=FluencebolHypothesis(HealpixNeutrinoAlertCatalogue()),
        tw=True
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_data_weight"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_data_weight.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "IR Flares",
        f"{len(catalogue_obj.parse_data())} sources",
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_data"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_data.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "IR Flares",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_scramble_1"),
        dpi=150,
        bbox_inches="tight"
    )
    plot_catalogue(
        catalogue_obj,
        5,
        "IR Flares",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_scramble_2"),
        dpi=150,
        bbox_inches="tight"
    )
    plot_catalogue(
        catalogue_obj,
        5,
        "IR Flares",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_scramble_3"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Nuclear Flares",
        "Flaires Catalog",
        nside=16,
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_binning"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Nuclear Flares",
        "Flaires Catalog",
        final_pdf=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"flaires_pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()

    # Fermi catalogue
    catalogue_obj = Fermi4FGLBlazarCatalogue()
    plot_catalogue(
        catalogue_obj,
        5,
        "Fermi Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        weighter=AverageFluxWeightHypothesis(
            HealpixNeutrinoAlertCatalogue()
        ),
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_data_weight"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_data_weight.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Fermi Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_data"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_data.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Fermi Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_scramble_1"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Fermi Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_scramble_2"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Fermi Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"fermi_scramble_3"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()

    from agn_radio_catalogue import AverageFluxWeightHypothesis

    # RFC catalogue
    catalogue_obj = AstrogeoAGNCatalogue()
    plot_catalogue(
        catalogue_obj,
        5,
        "RFC Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        weighter=AverageFluxWeightHypothesis(
            HealpixNeutrinoAlertCatalogue()
        ),
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_data_weight"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_data_weight.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "RFC Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_data"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_data.pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "RFC Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_scramble_1"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "RFC Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_scramble_2"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "RFC Blazars",
        f"{len(catalogue_obj.parse_data())} sources",
        scramble=True,
        scramble_size=8,
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_scramble_3"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Radio blazars",
        "Radio Catalog",
        nside=16,
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_binning"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
    plot_catalogue(
        catalogue_obj,
        5,
        "Radio blazars",
        "Radio Catalog",
        final_pdf=True,
    )
    plt.savefig(
        os.path.join(figures_folder,"rfc_pdf"),
        dpi=150,
        bbox_inches="tight"
    )
    plt.close()
