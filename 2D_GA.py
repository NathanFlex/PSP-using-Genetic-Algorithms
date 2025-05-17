import random
from typing import List, Tuple
from itertools import combinations
import tkinter as tk
import matplotlib.pyplot as plt
import pandas as pd

def measure_stability_directions(directions, sequence):
    # Convert to coordinates first
    coord_sequence = convert_to_coordinates(directions, sequence)
    if coord_sequence is None:  # Invalid structure
        return -1

    return len(measure_stability_coordinates(coord_sequence))

def measure_stability_coordinates(sequence):
    coords = [item[0] for item in sequence]
    labels = [item[1] for item in sequence]

    indices = [i for i, label in enumerate(labels) if label == 'H']
    points = [coords[i] for i in indices]

    adjacent_Hs = []
    for (idx1, p1), (idx2, p2) in combinations(zip(indices, points), 2):
        # Check if points are adjacent in grid but not in sequence
        if ((abs(p1[0] - p2[0]) == 1 and p1[1] == p2[1]) or
                (abs(p1[1] - p2[1]) == 1 and p1[0] == p2[0])):
            # Make sure they're not sequential in the chain
            if (p1, p2) not in zip(coords, coords[1:]) and (p2, p1) not in zip(coords, coords[1:]):
                adjacent_Hs.append((p1, p2))

    return adjacent_Hs
def convert_to_coordinates(directions, sequence):
    moves = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}
    coords = [(0, 0)]
    result = [((0, 0), sequence[0])]

    x, y = 0, 0
    for i, direction in enumerate(directions):
        dx, dy = moves[direction]
        x, y = x + dx, y + dy
        if (x, y) in coords:
            return None
        coords.append((x, y))
        result.append(((x, y), sequence[i + 1]))

    return result


def initial_population(sequence, pop_size = 20):
    population = []
    n = len(sequence) - 1

    while len(population) < pop_size:
        # directions = [random.randint(0, 3) for _ in range(n)]
        directions = [0 for _ in range(n)]

        if convert_to_coordinates(directions, sequence) is not None:
            population.append(directions)

    return population


def crossover(parent1, parent2, sequence):
    for _ in range(10):  # Try up to 10 times
        # Single point crossover
        point = random.randint(1, len(parent1) - 1)
        child = parent1[:point] + parent2[point:]

        # Is child valid?
        if convert_to_coordinates(child, sequence) is not None:
            return child

    # Return parents if not
    return random.choice([parent1, parent2]).copy()


def mutate(individual, sequence, mutation_rate = 0.2):
    for i in range(len(individual)):
        if random.random() < mutation_rate:
            original_direction = individual[i]
            possible_directions = [d for d in range(4) if d != original_direction]

            # Try all directions
            random.shuffle(possible_directions)
            for new_direction in possible_directions:
                mutated = individual.copy()
                mutated[i] = new_direction
                if convert_to_coordinates(mutated, sequence) is not None:
                    return mutated

    return individual


def evolve_population(population, sequence, select_style, select_size = 2):
    fitness_scores = [measure_stability_directions(ind, sequence) for ind in population]
    # print(population)
    # print(fitness_scores)


    if select_style != "Tournament":

        sorted_pop = [x for _, x in sorted(zip(fitness_scores, population),
                                           key=lambda p: p[0], reverse=True)]
        # print(sorted_pop)
        # print(sorted_pop[0:select_size])

        # Keep good individuals
        new_population = sorted_pop[0:select_size]
        # print(new_population)

        parent1, parent2 = new_population[0],new_population[1]

        while len(new_population) < len(population):
            # Crossover and mutation
            child = crossover(parent1, parent2, sequence)
            child = mutate(child, sequence)
            new_population.append(child)

    # Tournament style selection
    else:
        new_population = []
        # Make new population
        while len(new_population) < len(population):
            tournament_size = 5
            parent1 = max(random.sample(list(enumerate(population)), tournament_size),
                          key=lambda x: fitness_scores[x[0]])[1]
            parent2 = max(random.sample(list(enumerate(population)), tournament_size),
                          key=lambda x: fitness_scores[x[0]])[1]

            # Crossover and mutation
            child = crossover(parent1, parent2, sequence)
            child = mutate(child, sequence)
            new_population.append(child)

    return new_population
def draw_structure(coords):

    Tk = tk.Tk()
    Tk.title("Protein Folding Visualization")
    width = 1200
    height = 600
    padding = 50
    point_radius = 7
    Canvas = tk.Canvas(Tk, width=width, height=height, bg='white')
    Canvas.pack(expand=True, fill='both')
    # Calculate scaling
    x_coords = [p[0][0] for p in coords]
    y_coords = [p[0][1] for p in coords]
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    scale_x = (width - 2 * padding - 20) / max(1, max_x - min_x)
    scale_y = (height - 2 * padding - 50) / max(1, max_y - min_y)
    scale = min(scale_x, scale_y)

    center_x = width // 2
    center_y = (height - 50) // 2

    points = []

    # Draw H-H contacts
    contacts = measure_stability_coordinates(coords)
    for (x1, y1), (x2, y2) in contacts:
        cx1 = center_x + (x1 - (max_x + min_x) / 2) * scale
        cy1 = center_y + (y1 - (max_y + min_y) / 2) * scale
        cx2 = center_x + (x2 - (max_x + min_x) / 2) * scale
        cy2 = center_y + (y2 - (max_y + min_y) / 2) * scale

        Canvas.create_line(cx1, cy1, cx2, cy2,
                                fill='red', dash=(2, 2))

    #Draw points
    for i, ((x, y), acid) in enumerate(coords):
        canvas_x = center_x + (x - (max_x + min_x) / 2) * scale
        canvas_y = center_y + (y - (max_y + min_y) / 2) * scale
        points.append((canvas_x, canvas_y))

        color = 'gold' if acid == 'H' else 'blue'
        Canvas.create_oval(
            canvas_x - point_radius, canvas_y - point_radius,
            canvas_x + point_radius, canvas_y + point_radius,
            fill=color, outline='black'
        )

        if i > 0:
            Canvas.create_line(
                points[i - 1][0], points[i - 1][1],
                canvas_x, canvas_y,
                fill='black', width=2
            )

    # Show energy
    energy = -len(contacts)
    Canvas.create_text(
        20, 70,
        text=f"Energy = {energy}",
        anchor='w'
    )

df = pd.read_excel('sequences from Uniref50 database.xlsx')

print("Protein Structure Prediction Algorithm")
choice = int(input("Enter Protein Choice: "))

result = df[df['Number'] == choice]
# print(result)

values = result.iloc[0]
print("Chosen Protein AC/ID:",values['Protein AC/ID'])
print("Sequence Length:",values['Length'])
print("Known Lowest Energy:",values['Known Lowest Energy'])
print()

sequence = values['HP Sequence']
fitness_limit = -values['Known Lowest Energy']


# Initialize population
population = initial_population(sequence)
# print(population)
all_population = []
best_performers = []

# Track best solution
best_fitness = -1
best_solution = None

# Evolve for some generations
generation = 0
gen_counter = 1
current_best = 0

# Parameters
selection_style = "Tournament"
gen_limit = 10000

# 1. Algorithm Reset
gen_reset_enabled = False

# 2. Previous MSC
prev_msc = False
mut_rate = 0.2


if prev_msc == False:
    NC_count = 0  # No Change count
    while best_fitness < fitness_limit:
        if generation > gen_limit:
            print("Max number of generations reached")
            break
        min_fitness = 1000
        max_fitness = -1000
        population = evolve_population(population, sequence, selection_style, 2)

        # Check current best
        current_best = max(population,
                           key=lambda ind: measure_stability_directions(ind, sequence))
        current_fitness = measure_stability_directions(current_best, sequence)

        # Update Counters
        generation += 1
        gen_counter += 1
        NC_count += 1

        if current_fitness > best_fitness:
            best_fitness = current_fitness
            best_solution = current_best
            print(f"Generation {generation}: New best fitness = {-best_fitness}")
            gen_counter = 0             # If improvement, then reset counter
            NC_count = 0

        best_performers.append((best_solution,best_fitness))
        all_population.append(population)

        if gen_reset_enabled == True:
            if gen_counter == 1000:
                print("No change for 1000 generations, retrying")
                gen_counter = 0
                best_fitness = 0
                population = initial_population(sequence)

        if NC_count == 1000:
            print("Mutation rate updated")
            NC_count = 0
            mut_rate += 0.05


else:
    best_prev = [generation,population]
    NC_count = 0        # No Change count
    while best_fitness < fitness_limit:
        if generation > gen_limit:
            print("Max number of generations reached")
            break

        min_fitness = 1000
        max_fitness = -1000
        population = evolve_population(population, sequence, selection_style, 2)

        # Check current best
        current_best = max(population,
                           key=lambda ind: measure_stability_directions(ind, sequence))
        current_fitness = measure_stability_directions(current_best, sequence)

        if current_fitness > best_fitness:
            # print(population)
            best_fitness = current_fitness
            best_solution = current_best
            print(f"Generation {generation}: New best fitness = {-best_fitness}")
            best_prev = [generation,population]    # Assume newly found population as the best previous one
            gen_counter = 0     # If improvement, then reset counter
            # mut_rate = 0.2

        best_performers.append((best_solution, best_fitness))
        all_population.append(population)
        generation += 1
        gen_counter += 1

        if gen_counter == 1000:
            print("No change for 1000 generations, going back to previous best, generation - ",best_prev[0])
            NC_count += 1
            # print(best_prev)
            gen_counter = 0
            best_prev_mutated = []
            for i in best_prev[1]:          # Mutates each member of the previous population
                best_prev_mutated.append(mutate(i, sequence, mut_rate))
            population = best_prev_mutated
            # print(population)

        if NC_count == 1000:
            print("Mutation rate updated")
            NC_count = 0
            mut_rate += 0.1


print("Number of generations taken: ", generation)
best_coords = convert_to_coordinates(best_solution, sequence)
draw_structure(best_coords)
print(f"Final solution has {best_fitness} H-H contacts")
# print(best_performers)
min_fitness_values = [min([measure_stability_directions(ind, sequence) for ind in population]) for population in all_population]
max_fitness_values = [max([measure_stability_directions(ind,sequence) for ind in population]) for population in all_population]

# print(len(min_fitness_values))
# print(len(max_fitness_values))
generations_list = range(1, len(best_performers) + 1)
best_fitness_vals = [fit[1] for fit in best_performers]
# print(len(generations_list))
# print(len(best_performers))

fig, ax = plt.subplots()
ax.plot(generations_list, best_fitness_vals, label='Best Fitness', color='black')
ax.fill_between(generations_list, min_fitness_values, max_fitness_values, color='gray', alpha=0.3,label='Fitness Range')
ax.set_xlabel('Generation')
ax.set_ylabel('Fitness')
ax.set_title('Fitness Over Generations')
ax.legend()
plt.show()