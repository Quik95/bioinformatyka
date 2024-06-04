import math
import pickle
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Tuple

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

    def __add__(self, other: "HybridisationSequence"):
        return HybridisationSequence(self.sequence + other.sequence)


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
        source: str = "https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n=32&k=4&mode=alternating&intensity=0&position=0&sqpe=0&sqne=0&pose=0") -> HybridisationInstance:
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


def main():
    instance = load_from_file()
    print("Instance loaded, beginning calculations")

    even_start = instance.start[:]
    odd_start = instance.start[1:]

    graph = defaultdict(list)
    for i in range(0, len(instance.pairing_sequences)):
        for j in range(0, len(instance.pairing_sequences)):
            if i == j:
                continue
            overlappings = instance.pairing_sequences[i].find_all_matches(instance.pairing_sequences[j])
            if len(overlappings) > 0:
                for overlap in overlappings:
                    graph[instance.pairing_sequences[i]].append((instance.pairing_sequences[j], overlap))

    with open("/tmp/graph.gv", "w") as f:
        f.write("digraph {\n")
        for key, value in graph.items():
            for v in value:
                f.write(f'{key.sequence} -> {v[0].sequence} [label={v[1]}];\n')
        f.write("}")

    # even_result, odd_result = instance.get_starting_sequences()
    # instance.pairing_sequences.remove(odd_result)
    #
    # suffix_length = len(
    #     instance.pairing_sequences[0].sequence) - 2
    #
    # while len(instance.pairing_sequences) > 0:
    #     even_suffix = even_result[-suffix_length:]
    #     odd_suffix = odd_result[-suffix_length:]
    #
    #     even_match = instance.find_best_pairing(even_suffix, odd_suffix)
    #     instance.pairing_sequences.remove(even_match.sequence)
    #
    #     odd_match = instance.find_best_pairing(odd_suffix, even_suffix)
    #     instance.pairing_sequences.remove(odd_match.sequence)
    #
    #     assert even_match.score == odd_match.score
    #
    #     even_result = even_result + even_match.sequence[-2:]
    #     odd_result = odd_result + odd_match.sequence[-2:]
    #
    # assert len(even_result.sequence) == len(odd_result.sequence)
    # result = ""
    # for i in range(0, len(even_result.sequence), 2):
    #     result += even_result.sequence[i] + odd_result.sequence[i]
    #
    # print(result)


if __name__ == "__main__":
    main()
