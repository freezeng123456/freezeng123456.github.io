# Freezeng Knowledge Base Authoring Rules

These rules apply to all content and visual work in this repository.

## SVG-First Diagram Policy

1. All explanatory diagrams must be authored as SVG by default. This includes flowcharts, architecture diagrams, sequence diagrams, state machines, decision trees, memory-layout diagrams, algorithm overviews, and conceptual illustrations.
2. Do not introduce Mermaid diagrams, screenshots of diagrams, or rasterized text diagrams. Raster images remain appropriate only when the source is inherently raster-based, such as a photograph, UI screenshot, or experimental heatmap.
3. Draw each SVG as a deliberately designed infographic rather than a direct rendering of boxes and arrows. Use visual hierarchy, grouping, lanes, cards, restrained accent colors, and concise labels to communicate the underlying idea.
4. Keep diagram sources maintainable. Reuse and extend `scripts/generate-diagram-svgs.mjs` when a diagram belongs to the existing visual system.
5. Store bilingual diagram pairs under:
   - `content/assets/diagrams/<topic>/zh/<name>.svg`
   - `content/assets/diagrams/<topic>/en/<name>.svg`
6. Keep the Chinese and English versions structurally equivalent while rewriting labels naturally for each language.
7. Embed SVGs with standard Markdown image syntax and a content-root asset path:

   ```md
   ![Descriptive alternative text](assets/diagrams/<topic>/<lang>/<name>.svg)
   ```

   Do not use Quartz wiki embeds for SVG files because they may render as `<object>` elements instead of responsive images.

8. Every SVG must:
   - define a `viewBox`;
   - include an accessible `<title>` and `<desc>`;
   - use `role="img"` with the appropriate accessible references;
   - avoid embedded scripts, remote resources, and external font dependencies;
   - remain readable at the article content width without horizontal overflow.
9. Use a lightweight connector system: open chevron arrowheads, rounded caps and joins, restrained line weights, and arrow colors that match their paths. Avoid oversized filled-triangle arrowheads and visually heavy connectors.

## Adapting KM Articles and Figures

1. KM articles and figures may be used as reference material when the user has authorized access to them.
2. Extract the source's claims, reasoning, system relationships, visual metaphor, and information hierarchy, then reorganize them into an independent explanation.
3. Paraphrase prose and captions. Do not copy extended wording verbatim.
4. Figures may closely follow the source's conceptual structure when that structure is necessary to explain the idea, but every public SVG must be redrawn from scratch with original geometry, typography, wording, and composition. Do not trace, screenshot, or reproduce a source figure pixel for pixel.
5. Preserve technical meaning and important parameters. Do not invent values, results, causal relationships, or implementation details that are not supported by the source.
6. When attribution is appropriate and safe to publish, cite the source article or the underlying public reference near the adapted discussion.
7. This repository is publicly deployed. Before adapting internal material, exclude anything that is confidential, access-controlled, security-sensitive, or not authorized for public release. Never publish credentials, internal IP addresses, private endpoints, personal data, unreleased business metrics, proprietary source code, or operational security details.

## Diagram Quality Gate

Before publishing a diagram change:

1. Confirm that no Mermaid blocks remain in the affected pages.
2. Validate every generated SVG with `xmllint --noout`.
3. Run `npm run check`, `npm test`, and `npx quartz build`.
4. Confirm that Quartz renders SVGs as responsive `<img>` elements, not `<object>` elements.
5. Check both language versions in a real browser for loading failures, clipped labels, overlap, and horizontal overflow.
6. After deployment, verify that every affected page and SVG asset returns HTTP 200.

## Cursor Cloud specific instructions

This repo is a Quartz v5 static-site generator (no backend/database). Standard commands live in `README.md` and `package.json` scripts; use those. Notes below are only the non-obvious gotchas.

- Node is pinned to `v22.16.0` (`.node-version`) and provided via nvm's `default` alias; login shells pick it up automatically. If a non-login shell resolves a different `node`, run `source ~/.nvm/nvm.sh && nvm use default` first.
- `npx quartz plugin install` (and the `prebuild`/`install-plugins` hook) is a no-op here: it prints `⚠ No quartz.lock.json found` and is NOT required. All plugins are regular npm `@quartz-community/*` dependencies, so `npm ci` is sufficient and `npx quartz build` works directly.
- Run/preview the site with `npx quartz build --serve` (HTTP on port `8080`, WebSocket hot-reload on `3001`). `--serve` auto-enables `--watch`, so edits under `content/` trigger an incremental rebuild automatically — no restart needed.
- A brand-new `content/*.md` file logs `isn't yet tracked by git, dates will be inaccurate` until committed; this is expected and harmless.
