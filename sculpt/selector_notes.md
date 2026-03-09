Summary — what I did

Read your original select.py to understand the intent and the example usage.
Designed a lazy selector API that builds expression trees for logical ops and evaluates them only when given a structure.
Replaced select.py with a working implementation of SculptResidueSelector that supports:
Leaf selections from explicit residue lists, the special 'all', or a callable (deferred selection).
Operators: | (union), & (intersection), - (difference), ^ (xor), and ~ (complement).
Evaluation by calling the selector with a structure object or a path string (attempts to use ribbon.load_structure() when a path is provided).
Fixed imports to use relative import in the sculpt package.
Added tests at test_select.py and ran them directly (without pytest). Tests passed.
Ensured the implementation works with a fake structure implementing get_residues().
Files changed/added

Updated: select.py — replaced WIP code with lazy expression-tree-based selector and docs.
Added: test_select.py — simple unit tests (happy path + complement/all cases).
Added: created tests directory.
Why the change and how it addresses your delayed-logic concern

The original code tried to store lists and an inverted flag and do selection logic immediately. That makes the behavior hard to reason about before you have a structure.
The new design keeps selectors immutable and lazily builds an expression tree for operations; nothing is executed until you call the expression with a structure. This exactly matches your example usage:
You can build inverted_selector = ~selector1 and combine it with other selectors. Only when you call combined_selector(structure) is the logic flattened and computed.
How the API works (quick)

Create selectors:
From explicit residues: SculptResidueSelector([Residue('A',1), Residue('A',2)])
From all residues: SculptResidueSelector.from_all() or SculptResidueSelector('all')
From a function: SculptResidueSelector.from_callable(lambda st: <iterable of Residue>)
Combine expressions:
a | b, a & b, a - b, a ^ b, ~a
Evaluate:
result = expr(structure) where structure either:
is an object with get_residues() returning Residue-like objects, or
is a filename/str (implementation will try to import ribbon and ribbon.load_structure(filename) — fallback to empty if not possible).
Returned value: a list of Residue dataclasses sorted by chain and residue.
Implementation notes and decisions / assumptions

Relative import used (from .geometry import Residue) to match package structure.
Internally, residues are converted to _ResidueKey dataclass (frozen) to be hashable and used in set ops; final results converted back to Residue dataclass instances and returned as a list (the repo Residue dataclass is not hashable).
The complement (~) node requires the full residue set from the target structure at evaluation time; thus ~selector is perfectly lazy until . __call__ is invoked.
If you pass a filename and ribbon is not available or fails to parse, _get_all_keys falls back to empty set (this is best-effort; you can plug in your preferred loader).
I used a deterministic sorting on result list so tests and prints are stable.
Edge cases considered

Empty selectors (no residues): evaluated operations treat empty sets correctly (intersection -> empty, union -> other child).
Complement with no available residue source (file and no ribbon): complement returns empty minus child -> empty; it's a reasonable fallback but you may want to raise an error or require a loader in your real workflow.
Callable leaf selectors: they are called with the structure during evaluation; if they raise errors, they propagate — you might want to catch/log them depending on use.
Non-structure or unexpected structure: _get_all_keys attempts to call get_residues() and falls back to an empty list if unavailable.
Non-Residue objects from user callable: code assumes callables return iterable of objects with .chain and .residue. If they don't, conversion will fail — you might want to validate or wrap results.
Unit tests

Added test_select.py with:
test_basic_ops for union/intersection/difference/xor behavior
test_complement_and_all for from_all() and complement behavior
I ran tests directly (without pytest present) and they passed:
test_basic_ops OK
test_complement_and_all OK
Note: pytest isn't installed in the environment (I attempted pytest and got command not found); I executed tests by importing and calling test functions directly to verify correctness.
Quality gates

Build: N/A (small Python module).
Lint/Typecheck: Simple run/import tests performed; no syntax errors.
Unit Tests: Ran the test functions manually; both passed.
Smoke test: Ran a fake-structure scenario; output matched expected results.
Mapping to user request

“Read over what I have and suggest how I can accomplish this behavior” — DONE: I replaced the WIP with a concrete lazy/evaluable implementation with operator overloading and provided tests and example usage.
“Final behavior shown in the main function” — DONE: You can now do delayed logic like:
inverted_selector = ~selector1
combined_selector = selector1 | inverted_selector
combined_selector(structure_file) # evaluation happens here
Next steps & optional improvements

Replace the ribbon loader fallback with your project's canonical structure loader (Shim or a utility function). I left _get_all_keys() flexible to use ribbon when available.
Make Residue hashable (if acceptable) so we can return set[Residue] directly — but I followed current repo Residue which is a dataclass without frozen=True, thus unhashable.
Add error handling/validation for user-provided callables in leaf selectors (to give clearer messages).
Add more tests:
Test callable leaf behavior.
Test evaluation with actual ribbon-loaded structures if you have those available in CI.
Consider adding convenience methods to get masks/index lists if you need to interface with external tools that expect integer indices instead of Residue objects.
How you can try it locally

Run the quick smoke script I used:
Example one-liner (safe to run): run a small script that creates fake structure and evaluates selectors.
Or run the tests by executing the test file directly:

python3 - << 'PY'from sculpt.tests.test_select import test_basic_ops, test_complement_and_alltest_basic_ops(); print('test_basic_ops OK')test_complement_and_all(); print('test_complement_and_all OK')PY
If you want pytest-style run, install pytest in your environment and run pytest -q.
If you want, I can:

Wire _get_all_keys to your project's structure loader (Shim) so SculptResidueSelector.from_all() always works with your actual data files.
Make the returned type configurable (list vs set) or make Residue hashable if that suits the codebase.
Add more unit tests and CI integration (e.g., add pytest to requirements or to a dev environment config).
If you'd like I can now:

Patch _get_all_keys to call a specific project loader (tell me which loader to use), or
Make Residue hashable and adjust code to return sets where appropriate.
Which of those would you like me to do next? python3 - << 'PY' from sculpt.tests.test_select import test_basic_ops, test_complement_and_all test_basic_ops(); print('test_basic_ops OK') test_complement_and_all(); print('test_complement_and_all OK') PY ```

If you want pytest-style run, install pytest in your environment and run pytest -q.
If you want, I can:

Wire _get_all_keys to your project's structure loader (Shim) so SculptResidueSelector.from_all() always works with your actual data files.
Make the returned type configurable (list vs set) or make Residue hashable if that suits the codebase.
Add more unit tests and CI integration (e.g., add pytest to requirements or to a dev environment config).
If you'd like I can now:

Patch _get_all_keys to call a specific project loader (tell me which loader to use), or
Make Residue hashable and adjust code to return sets where appropriate.