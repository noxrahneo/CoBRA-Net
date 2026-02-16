# Git Hooks

This folder contains shared git hook scripts for the repository.

## Enable Hooks Locally

Git does not automatically load hooks from versioned folders. To enable these hooks on my machine:

```bash
git config core.hooksPath hooks
```

## Available Hooks

- `pre-commit`: Updates the **Last updated** date in the top-level README.md.

## How It Works

1. Run `git config core.hooksPath hooks` once per machine/repo.
2. On every `git commit`, Git runs `hooks/pre-commit` automatically.
3. The hook updates the **Last updated** line in README.md and stages it.

## Notes

- Hooks are executed locally. Teammates must run the configuration command once.
- You can add more hook scripts here as needed (e.g., `pre-push`, `commit-msg`).
