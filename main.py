import heapq
import math
import pickle
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Tuple, Set, Union

import requests
from lxml import etree


@dataclass(frozen=True, eq=True)
class HybridisationSequence:
    sequence: str

    def __getitem__(self, item):
        return HybridisationSequence(self.sequence[item])

    def __add__(self, other: Union["HybridisationSequence", str]):
        if isinstance(other, str):
            return HybridisationSequence(self.sequence + other)

        return HybridisationSequence(self.sequence + other.sequence)

    def __len__(self) -> int:
        return len(self.sequence)

    def _overlap(self, other: "HybridisationSequence") -> bool:
        assert self[0] != "X" and other[0] != "X"
        return [c for c in self.sequence if c != "X"] == [c for c in other.sequence[:len(self)] if c != "X"]

    def get_overlaps(self, other: "HybridisationSequence") -> List[int]:
        overlaps = []
        for i in range(0, len(self), 2):
            if self[i:]._overlap(other):
                overlaps.append(i // 2)

        # Due to malformed input, the even starting vertex always is smaller than the rest by two characters
        # if that is the case, we add one two each overlap to account for that
        if len(other) - len(self) == 2:
            overlaps = [o + 1 for o in overlaps]

        return overlaps


@dataclass(frozen=True, eq=True)
class GraphArc:
    start: HybridisationSequence
    end: HybridisationSequence
    weight: int


@dataclass
class HybridisationInstance:
    start: HybridisationSequence
    target_length: int
    pairing_sequences: List[HybridisationSequence]
    testing_sequences: Set[HybridisationSequence]
    l: int

    @staticmethod
    def load_from_generator(n: int = 32, k: int = 4, sqnep: int = 0) -> "HybridisationInstance":
        url = f"https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n={n}&k={k}&mode=alternating&intensity=0&position=0&sqpe=0&sqnep={sqnep}&pose=0"
        r = requests.get(url)
        assert r.status_code == 200
        if not r.text.startswith("<"):
            raise ValueError("Invalid response", r.text)

        return HybridisationInstance._load_instance(r.text)

    @staticmethod
    def load_from_file(source: str = "./small.xml") -> "HybridisationInstance":
        with open(source, "r") as f:
            return HybridisationInstance._load_instance(f.read())

    @staticmethod
    def _load_instance(source: str) -> "HybridisationInstance":
        root = etree.fromstring(source)

        start = HybridisationSequence(root.get("start"))
        target_length = int(root.get("length"))
        pairing_sequences = [HybridisationSequence(el.text) for el in root.getchildren()[0]]
        testing_sequences = set([HybridisationSequence(el.text) for el in root.getchildren()[1]])

        return HybridisationInstance(
            start=start,
            target_length=target_length,
            pairing_sequences=pairing_sequences,
            testing_sequences=testing_sequences,
            l=len(start)
        )

    def get_staring_vertices(self) -> Tuple[HybridisationSequence, HybridisationSequence]:
        return (
            HybridisationSequence("".join([self.start.sequence[i] + "X" for i in range(0, len(self.start), 2)])[:-1]),
            HybridisationSequence("".join([self.start.sequence[i] + "X" for i in range(1, len(self.start), 2)]))[:-1])

    def verify_vertex(self, candidate: HybridisationSequence, complementary: HybridisationSequence) -> bool:
        return complementary[2:len(complementary) - 1] + candidate[-1] in self.testing_sequences

    def is_solution_complete(self, solution: "PartialSolution") -> bool:
        return len(solution.odd_sequence) + 1 == self.target_length and len(
            solution.even_sequence) + 1 == self.target_length

    def can_be_pruned(self, solution: "PartialSolution") -> bool:
        return len(solution.odd_sequence) + 1 > self.target_length or len(
            solution.even_sequence) + 1 > self.target_length


class PartialSolution:
    visited_vertices: Set[HybridisationSequence]
    odd_sequence: HybridisationSequence
    even_sequence: HybridisationSequence

    def __init__(self, visited_vertices: Set[HybridisationSequence], odd_sequence: HybridisationSequence,
                 even_sequence: HybridisationSequence):
        self.visited_vertices = visited_vertices
        self.odd_sequence = odd_sequence
        self.even_sequence = even_sequence

    def extend(self, odd_arc: GraphArc, even_arc: GraphArc) -> "PartialSolution":
        new_visited = self.visited_vertices.union({odd_arc.end, even_arc.end})
        new_odd = self.odd_sequence + odd_arc.end[-odd_arc.weight * 2:]
        new_even = self.even_sequence + even_arc.end[-even_arc.weight * 2:]

        return PartialSolution(new_visited, new_odd, new_even)

    def current_vertices(self, l: int) -> Tuple[HybridisationSequence, HybridisationSequence]:
        return self.odd_sequence[-l:], self.even_sequence[-max(0, l):]

    def reconstruct(self) -> str:
        result = ""
        for i in range(0, len(self.odd_sequence), 2):
            result += self.odd_sequence.sequence[i] + self.even_sequence.sequence[i]
        return result


class HeapItem:
    def __init__(self, item: PartialSolution):
        self.priority = -len(item.visited_vertices)
        self.item = item

    def value(self) -> PartialSolution:
        return self.item

    def __lt__(self, other):
        return self.priority < other.priority

    def __eq__(self, other):
        return self.priority == other.priority


def main():
    # instance = HybridisationInstance.load_from_file()
    instance = HybridisationInstance.load_from_generator(n=1000, k=10, sqnep=1)
    print("Instance loaded, beginning calculations")

    odd_start, even_start = instance.get_staring_vertices()
    graph = defaultdict(list)

    p = instance.pairing_sequences + [odd_start,
                                      even_start] if odd_start not in instance.pairing_sequences else instance.pairing_sequences + [
        even_start]
    for start in p:
        for end in instance.pairing_sequences:
            if start != end:
                for overlap in start.get_overlaps(end):
                    graph[start].append(GraphArc(start, end, overlap))

    for start, arcs in graph.items():
        graph[start] = sorted(arcs, key=lambda x: x.weight)

    with open("/tmp/graph.gv", "w") as f:
        f.write("digraph G {\n")
        for start, arcs in graph.items():
            start_label = f"odd start\n{start.sequence}" if start == odd_start else start.sequence
            start_label = f"even start\n{start.sequence}" if start == even_start else start_label
            for arc in arcs:
                f.write(
                    f'"{start_label}" -> "{arc.end.sequence}" [label="{arc.weight}"] [color="{"red" if arc.weight == 1 else "black"}"];\n')
        f.write("}")

    print("Graph created, beginning calculations")
    initial = PartialSolution({odd_start, even_start}, odd_sequence=odd_start, even_sequence=even_start)

    queue = [HeapItem(initial)]
    heapq.heapify(queue)

    solutions = []
    timelimit = time.time() + 60 * 2  # 2 minutes
    while len(queue) > 0:
        if time.time() > timelimit:
            print("Time limit exceeded")
            break

        current = heapq.heappop(queue).value()
        current_odd_vertex, current_even_vertex = current.current_vertices(len(instance.start))

        odd_candidates = [arc for arc in graph[current_odd_vertex] if arc.end not in current.visited_vertices]
        even_candidates = [arc for arc in graph[current_even_vertex] if
                           arc.end not in current.visited_vertices and arc.end not in odd_candidates]

        for odd_arc in odd_candidates:
            if (odd_arc.weight == 1
                    and current_even_vertex != even_start  # Due to malformed input, in the first iteration the correct vertex would be rejected
                    and not instance.verify_vertex(odd_arc.end, current_even_vertex)):
                continue

            for even_arc in even_candidates:
                if even_arc.weight == 1 and not instance.verify_vertex(even_arc.end, current_odd_vertex):
                    continue

                new_solution = current.extend(odd_arc, even_arc)
                if instance.is_solution_complete(new_solution):
                    solutions.append(new_solution)
                elif not instance.can_be_pruned(new_solution):
                    heapq.heappush(queue, HeapItem(new_solution))

                # If we got here and the arc has weight 1, that means that it was verified and there can't any better solutions
                if even_arc.weight == 1:
                    break

            # If we got here and the arc has weight 1, that means that it was verified and there can't any better solutions
            if odd_arc.weight == 1:
                break

    print(f"Calculations finished. Number of unique solutions found: {len(set([s.reconstruct() for s in solutions]))}")
    pass


if __name__ == "__main__":
    main()
