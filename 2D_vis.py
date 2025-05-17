import tkinter as tk
import random
from typing import List, Tuple
from itertools import combinations
import time


class ProteinFolder:
    def __init__(self, sequence: str, population_size: int = 10, width=1200, height=600):
        self.root = tk.Tk()
        self.root.title("Protein Folding Visualization")

        # Canvas setup
        self.canvas = tk.Canvas(self.root, width=width, height=height, bg='white')
        self.canvas.pack(expand=True, fill='both')

        # Drawing parameters
        self.width = width
        self.height = height
        self.padding = 50
        self.point_radius = 5

        # Create subplot areas
        self.left_subplot = self.canvas.create_rectangle(
            10, 10, width // 2 - 10, height - 10, outline='gray')
        self.right_subplot = self.canvas.create_rectangle(
            width // 2 + 10, 10, width - 10, height - 10, outline='gray')

        # Labels
        self.canvas.create_text(width // 4, 30, text="Best Solution", font=('Arial', 14))
        self.canvas.create_text(3 * width // 4, 30, text="Current Population", font=('Arial', 14))

        # Evolution parameters
        self.sequence = sequence
        self.population_size = population_size
        self.population = initial_population(self.sequence)
        self.best_solution = None
        self.best_fitness = 0
        self.generation = 0
        self.current_member = 0

        # Animation control
        self.delay = 100
        self.is_running = False

        # Start/Stop button
        self.button = tk.Button(self.root, text="Start Evolution", command=self.toggle_evolution)
        self.button.pack()

        # Generation label
        self.gen_label = tk.Label(self.root, text="Generation: 0")
        self.gen_label.pack()

        # Fast forward Button
        self.ff = tk.Button(self.root,text="Fast Forward",command=self.switch_speed)
        self.ff.pack()

    def switch_speed(self):
        if self.delay == 100:
            self.delay = 10
            self.ff.config(text='Slow down')
        else:
            self.delay = 100
            self.ff.config(text='Fast Forward')

    def toggle_evolution(self):
        if not self.is_running:
            self.is_running = True
            self.button.config(text="Stop Evolution")
            self.run_evolution()
        else:
            self.is_running = False
            self.button.config(text="Start Evolution")

    def draw_structure(self, coords, subplot='left', clear=True):
        # Calculate subplot boundaries
        if subplot == 'left':
            x_offset = 10
            width = self.width // 2 - 20
        else:
            x_offset = self.width // 2 + 10
            width = self.width // 2 - 20

        if subplot == 'left':
            tag = "best"
        else:
            tag = "current"

        if clear:
            self.canvas.create_rectangle(
                x_offset, 40, x_offset + width, self.height - 10,
                fill='white', outline='gray',tags="tag")

            if subplot == 'right':
                self.canvas.delete("current")
            else:
                self.canvas.delete("best")

        # Calculate scaling
        x_coords = [p[0][0] for p in coords]
        y_coords = [p[0][1] for p in coords]
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)

        scale_x = (width - 2 * self.padding) / (max(1, max_x - min_x))
        scale_y = (self.height - 2 * self.padding - 50) / (max(1, max_y - min_y))
        scale = min(scale_x, scale_y)

        center_x = x_offset + width // 2
        center_y = (self.height - 50) // 2

        points = []
        for i, ((x, y), acid) in enumerate(coords):
            canvas_x = center_x + (x - (max_x + min_x) / 2) * scale
            canvas_y = center_y + (y - (max_y + min_y) / 2) * scale
            points.append((canvas_x, canvas_y))

            # Draw amino acid
            color = 'gold' if acid == 'H' else 'black'

            self.canvas.create_oval(
                canvas_x - self.point_radius, canvas_y - self.point_radius,
                canvas_x + self.point_radius, canvas_y + self.point_radius,
                fill=color, outline='black', tags=tag)

            if i > 0:
                self.canvas.create_line(
                    points[i - 1][0], points[i - 1][1],
                    canvas_x, canvas_y,
                    fill='blue', width=2, tags=tag)

        # Draw H-H contacts
        contacts = measure_stability_coordinates(coords)
        for (x1, y1), (x2, y2) in contacts:
            # Transform coordinates
            cx1 = center_x + (x1 - (max_x + min_x) / 2) * scale
            cy1 = center_y + (y1 - (max_y + min_y) / 2) * scale
            cx2 = center_x + (x2 - (max_x + min_x) / 2) * scale
            cy2 = center_y + (y2 - (max_y + min_y) / 2) * scale

            self.canvas.create_line(cx1, cy1, cx2, cy2,
                                    fill='red', dash=(2, 2), tags = tag)

        # Show energy
        energy = -len(contacts)
        self.canvas.create_text(
            x_offset + 20, 70,
            text=f"Energy = {energy}",
            anchor='w', tags = tag)

    def run_evolution(self):
        if not self.is_running:
            return

        self.gen_label.config(text=f"Generation: {self.generation}")

        coords = convert_to_coordinates(
            self.population[self.current_member],
            self.sequence)
        if coords:
            self.draw_structure(coords, 'right')

            # Update best solution
            fitness = -len(measure_stability_coordinates(coords))
            # print(fitness)
            # print(self.best_fitness)
            if fitness < self.best_fitness:
                self.best_fitness = fitness
                self.best_solution = self.population[self.current_member]
                best_coords = coords
                self.draw_structure(best_coords, 'left')

        # Move to next member
        self.current_member += 1

        # If we've shown all members, evolve to next generation
        if self.current_member >= len(self.population):
            self.population = evolve_population(
                self.population, self.sequence)
            self.generation += 1
            self.current_member = 0

        # Schedule next update
        self.root.after(self.delay, self.run_evolution)

    def start(self):
        self.root.mainloop()

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


def convert_to_directions(sequence):
    directions = []
    for i in range(1, len(sequence)):
        curr_x, curr_y = sequence[i][0]
        prev_x, prev_y = sequence[i - 1][0]
        dx = curr_x - prev_x
        dy = curr_y - prev_y

        if dx == 1:
            directions.append(0)
        elif dx == -1:
            directions.append(2)
        elif dy == 1:
            directions.append(1)
        elif dy == -1:
            directions.append(3)

    return directions


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


def initial_population(sequence, pop_size = 10):
    population = []
    n = len(sequence) - 1

    while len(population) < pop_size:
        directions = [random.randint(0, 3) for _ in range(n)]

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


def evolve_population(population, sequence, select_size = 2):
    fitness_scores = [measure_stability_directions(ind, sequence) for ind in population]

    sorted_pop = [x for _, x in sorted(zip(fitness_scores, population),
                                       key=lambda p: p[0], reverse=True)]
    # print(sorted_pop)
    # print(sorted_pop[0:select_size])

    # Keep good individuals
    new_population = sorted_pop[0:select_size]
    # print(new_population)

    parent1, parent2 = new_population[0], new_population[1]

    while len(new_population) < len(population):
        # Crossover and mutation
        child = crossover(parent1, parent2, sequence)
        child = mutate(child, sequence)
        new_population.append(child)

    return new_population

protein = 'HHPHHHHHPPPPPH'
vis = ProteinFolder(protein)
vis.start()