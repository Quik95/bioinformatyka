import math
import pickle
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Tuple, Set, Union

import requests
from lxml import etree


@dataclass
class SortItem:
    sequence: "HybridisationSequence"
    score: int


@dataclass(frozen=True, eq=True)
class HybridisationSequence:
    sequence: str

    def find_all_matches(self, other: "HybridisationSequence") -> List[int]:
        """ Finds all indexes where two sequences begin to match """
        matches = []
        for i in range(0, len(self.sequence), 2):
            match = True
            for j, nucl in enumerate(self.sequence[i:]):
                if j >= len(other.sequence):
                    break
                if nucl == "X" or other.sequence[j] == "X":
                    continue
                if nucl != other.sequence[j]:
                    match = False
                    break
            if match:
                matches.append(i)
        return matches

    def find_match(self, other: "HybridisationSequence") -> Optional[int]:
        """ Finds the index where two sequences begin to match """
        for i in range(0, len(self.sequence), 2):
            match = True
            for j, nucl in enumerate(self.sequence[i:]):
                if j >= len(other.sequence):
                    break
                if nucl == "X" or other.sequence[j] == "X":
                    continue
                if nucl != other.sequence[j]:
                    match = False
                    break
            if match:
                return i
        return None

    @staticmethod
    def sort_matches(needle: "HybridisationSequence", haystack: List["HybridisationSequence"]) \
            -> List[SortItem]:
        return sorted(
            [SortItem(seq, HybridisationSequence(needle.sequence).find_match(seq)) for seq in haystack],
            key=lambda x: x.score if x.score is not None else math.inf
        )

    def __getitem__(self, item):
        return HybridisationSequence(self.sequence[item])

    def __add__(self, other: Union["HybridisationSequence", str]):
        if isinstance(other, str):
            return HybridisationSequence(self.sequence + other)

        return HybridisationSequence(self.sequence + other.sequence)


@dataclass(frozen=True, eq=True)
class GraphArc:
    start: str
    end: str
    overlap: int


@dataclass
class HybridisationInstance:
    start: HybridisationSequence
    target_length: int
    pairing_sequences: List[HybridisationSequence]
    testing_sequences: List[HybridisationSequence]

    def get_starting_sequences(self) -> Tuple[HybridisationSequence, HybridisationSequence]:
        even_starts = HybridisationSequence.sort_matches(self.start, self.pairing_sequences)
        assert even_starts[1].score != 0
        even_start = even_starts[0].sequence

        odd_starts = HybridisationSequence.sort_matches(self.start[1:], self.pairing_sequences)
        assert odd_starts[1].score != 0
        odd_start = odd_starts[0].sequence

        return even_start, odd_start

    def test_found_matches(self, matches: List[SortItem], other_sequence_suffix: HybridisationSequence) -> SortItem:
        for potential in matches:
            if other_sequence_suffix + potential.sequence[-1] in self.testing_sequences:
                return potential

    def verify_extension(self, extension: GraphArc, secondary: HybridisationSequence) -> bool:
        secondary_suffix = secondary[-len(extension.end) + extension.overlap:] + extension.end[-1]
        return secondary_suffix in self.testing_sequences

    def find_best_pairing(self, primary: HybridisationSequence, secondary: HybridisationSequence) -> SortItem:
        matches = [seq for seq in HybridisationSequence.sort_matches(primary, self.pairing_sequences) if
                   seq.score is not None]
        assert len(matches) > 0
        best_match = matches[0]
        if len([seq for seq in matches if seq.score == 0]) > 1:
            # there is more than one perfect match
            best_match = self.test_found_matches(matches, secondary)
        return best_match


def load_from_url(
        source: str = "https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n=16&k=4&mode=alternating&intensity=0&position=0&sqpe=0&sqne=0&pose=0") -> HybridisationInstance:
    r = requests.get(source)
    assert r.status_code == 200
    if not r.text.startswith("<"):
        raise ValueError("Invalid response", r.text)

    return load_instance(r.text)


def load_from_file(source: str = "./small.xml") -> HybridisationInstance:
    with open(source, "r") as f:
        return load_instance(f.read())


def load_instance(source: str) -> HybridisationInstance:
    root = etree.fromstring(source)

    start = HybridisationSequence(root.get("start"))
    target_length = int(root.get("length"))
    pairing_sequences = [HybridisationSequence(el.text) for el in root.getchildren()[0]]
    testing_sequences = [HybridisationSequence(el.text) for el in root.getchildren()[1]]

    return HybridisationInstance(
        start=start,
        target_length=target_length,
        pairing_sequences=pairing_sequences,
        testing_sequences=testing_sequences
    )


@dataclass
class PartialSolution:
    used_vertices: Set[str]
    current_even_vertex: str
    current_odd_vertex: str
    odd_sequence: HybridisationSequence
    even_sequence: HybridisationSequence

    @staticmethod
    def create_new(even_start: str, odd_start: str) -> "PartialSolution":
        return PartialSolution(
            used_vertices={even_start, odd_start},
            current_even_vertex=even_start,
            current_odd_vertex=odd_start,
            odd_sequence=HybridisationSequence(odd_start),
            even_sequence=HybridisationSequence(even_start),
        )

    def extend(self, even: GraphArc, odd: GraphArc) -> "PartialSolution":
        new_odd = HybridisationSequence(odd.end) if odd.overlap == 0 and len(odd.end) > len(
            self.odd_sequence.sequence) else self.odd_sequence + odd.end[-odd.overlap:]
        new_even = HybridisationSequence(even.end) if even.overlap == 0 and len(even.end) > len(
            self.even_sequence.sequence) else self.even_sequence + even.end[-even.overlap:]

        return PartialSolution(
            used_vertices=self.used_vertices.union({even.end, odd.end}),
            current_even_vertex=even.end,
            current_odd_vertex=odd.end,
            odd_sequence=new_odd,
            even_sequence=new_even,
        )

    def is_size_correct(self, target: int) -> bool:
        return len(self.odd_sequence.sequence) + 1 == target and len(self.even_sequence.sequence) + 1 == target

    def can_be_pruned(self, target: int) -> bool:
        return len(self.odd_sequence.sequence) + 1 > target or len(self.even_sequence.sequence) + 1 > target

    def recombine(self) -> str:
        result = ""
        for i in range(0, len(self.even_sequence.sequence), 2):
            result += self.even_sequence.sequence[i] + self.odd_sequence.sequence[i]
        return result


def main():
    instance = load_from_file()
    # instance = load_from_url()
    print("Instance loaded, beginning calculations")

    even_start = "".join([instance.start.sequence[i] + "X" for i in range(0, len(instance.start.sequence), 2)])[:-1]
    odd_start = "".join([instance.start.sequence[i] + "X" for i in range(1, len(instance.start.sequence), 2)])[:-1]
    instance.pairing_sequences += [HybridisationSequence(even_start), HybridisationSequence(odd_start)]

    graph = defaultdict(list)
    for i in range(0, len(instance.pairing_sequences)):
        for j in range(0, len(instance.pairing_sequences)):
            if i == j:
                continue
            overlappings = instance.pairing_sequences[i].find_all_matches(instance.pairing_sequences[j])
            if len(overlappings) > 0:
                for overlap in overlappings:
                    graph[instance.pairing_sequences[i].sequence].append(
                        GraphArc(start=instance.pairing_sequences[i].sequence,
                                 end=instance.pairing_sequences[j].sequence,
                                 overlap=overlap))

    with open("/tmp/graph.gv", "w") as f:
        f.write("digraph {\n")
        for key, value in graph.items():
            for v in value:
                f.write(f"{key} -> {v.end} [label={v.overlap}];\n")
        f.write("}")

    solutions = []
    processing_queue = [PartialSolution.create_new(even_start, odd_start)]
    while len(processing_queue) > 0:
        current = processing_queue.pop(0)
        if current.can_be_pruned(instance.target_length):
            continue

        if current.is_size_correct(instance.target_length):
            if len(current.odd_sequence.sequence) == len(current.even_sequence.sequence):
                solutions.append(current)
                continue

        even_candidates = [arc for arc in graph[current.current_even_vertex] if arc.end not in current.used_vertices]
        odd_candidates = [arc for arc in graph[current.current_odd_vertex] if
                          arc.end not in current.used_vertices and arc not in even_candidates]

        for even_candidate in even_candidates:
            break_loop = False
            # current verification is incorrect, rework it next
            if (even_candidate.overlap == 2
                    and current.current_odd_vertex[2:] + even_candidate.end[-1] in instance.testing_sequences):
                # if we can verify the shortest possible sequence, we can't add any other vertices as they are sure to
                # be invalid, however, if we can't verify the shortest possible sequence, can't discard it as it might
                # be a part of a valid sequence but its verifying sequence is missing due to presence of negative errors
                break_loop = True  # try to add the even candidate and then break the loop
            for odd_candidate in odd_candidates:
                new_solution = current.extend(even_candidate, odd_candidate)
                processing_queue.append(new_solution)

                if odd_candidate.overlap == 2:
                    if current.current_even_vertex[2:] + odd_candidate.end[-1] not in instance.testing_sequences:
                        # verification failed, continue
                        processing_queue.pop()
                        continue
                    else:
                        break

            if break_loop:
                break
    solutions[0].recombine()
    pass


if __name__ == "__main__":
    main()
