"""Compatibility entrypoint that delegates to bio_mcp.server."""

from bio_mcp.server import build_server, main

__all__ = ["build_server", "main"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
