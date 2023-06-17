fractions = {
    'Gas1': [
        [0.3, 0.2, 0.5],
        [0.1, 0.4, 0.5],
        [0.2, 0.3, 0.5]
    ],
    'Gas2': [
        [0.4, 0.5, 0.1],
        [0.2, 0.3, 0.5],
        [0.3, 0.1, 0.6]
    ]
}

first_key = list(fractions.keys())[0]
num_elements = len(fractions[first_key][0])

oil_fractions = []

for i in range(len(fractions[first_key])):
    for j in range(num_elements):
        oil_fraction = round(1 - sum(fractions[key][i][j] for key in fractions), 1)
        oil_fractions.append(oil_fraction)

print(oil_fractions)