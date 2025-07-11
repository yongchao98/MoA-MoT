# A list of League of Legends champions that can buffer their dash ability with Flash
# to redirect its primary damage or crowd-control effect.
champions_who_can_dash_flash = [
    "Gragas",
    "Vi",
    "Jarvan",
    "Galio",
    "Ornn",
    "Rakan",
    "Rell",
    "Riven",
    "Sejuani",
    "Shen",
    "Tryndamere",
    "Urgot",
    "Yone",
    "Zac",
    "Camille"
]

# The script will join the elements of the list into a single string,
# with each name separated by a comma.
output_string = ",".join(champions_who_can_dash_flash)

# Print the final comma-separated list of champion names.
print(output_string)