# MandrakeModels lessons

## Editing published files

- A line-delete (`sed '/pattern/d'`) to strip a private path also deletes whatever else shares that line; in
  `local_inertial_ripsense.py` it took the closing `"""` of a docstring and the module stopped importing
  (shipped in 75766e6, fixed in e17386b). Replace the fragment, then `python3 -m py_compile` before committing.
- This repo is PUBLIC (since 2026-11-15). New text names no private paths (`D:\ClaudeWork`), host names or tracker
  item ids; commits use the GitHub noreply address (repo-local `user.email`).

- 2026-11-16 (handoff SOLO_HANDOFF_202611161515): The repo is PUBLIC on github.com (MIT, since 2026-11-15): anything committed is published. Keep private paths, host names and tracker ids out of new files and commit messages.
