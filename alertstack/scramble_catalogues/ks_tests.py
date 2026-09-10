from agn_radio_catalogue import AstrogeoAGNCatalogue
from flaires_catalogue import FlairesCatalogue
from sjoert_catalogue import AccretionFlaresSjoertCatalogue


if __name__ == "__main__":

    print("\nKS Tests for the 63 Accretion Flares")
    catalogue_obj = AccretionFlaresSjoertCatalogue()
    catalogue_obj.mykstest_decs()
    catalogue_obj.mykstest_ras()

    print("\nKS Tests for the Flaires catalogue")
    catalogue_obj = FlairesCatalogue()
    catalogue_obj.mykstest_decs()
    catalogue_obj.mykstest_ras()

    print("\nKS Tests for the Radio Fundamental Catalogue")
    catalogue_obj = AstrogeoAGNCatalogue()
    catalogue_obj.mykstest_decs()
    catalogue_obj.mykstest_ras()