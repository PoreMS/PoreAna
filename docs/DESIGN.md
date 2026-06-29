# PoreAna Visual Design

## Color Palette

The logo and all visual assets use a strict 3-color red palette mirroring the
[Material Design Red scale](https://m2.material.io/design/color/the-color-system.html),
following the same structural logic as PoreMS (blue) and PoreSim (amber).

| Role | Hex | RGB | Material token |
|---|---|---|---|
| Bond lines | `#EF9A9A` | `rgb(239, 154, 154)` | Red 200 |
| Small/medium nodes, accent text ("ANA") | `#C74949` | `rgb(199, 73, 73)` | Red 500 |
| Large anchor nodes, body text ("PORE") | `#724242` | `rgb(114, 66, 66)` | Red 900 |

### Rationale

- **3 colors maximum** keeps the graphic legible at small sizes (favicon, sidebar logo).
- **Perceptual hierarchy** matches PoreMS exactly: bonds are the lightest element (connectors),
  large anchor nodes are the darkest (structural weight).
- **Material Red scale** provides perceptually uniform steps and is visually distinct from
  PoreMS (blue) and PoreSim (amber) while sharing the same design grammar.

## Logo Files

| File | Usage |
|---|---|
| `docs/pics/logo.svg` | Square icon (favicon source, app icon) |
| `docs/pics/logo_text.svg` | Horizontal logo with "PoreAna" wordmark |
| `docs/pics/logo_text_sub.svg` | Logo with wordmark + subtitle line |

## Typography

Logo text uses **Arial / Arial MT** (system sans-serif fallback).
- "**PORE**": Red 900 (`#724242`) — darkest, matches large anchor nodes
- "**ANA**": Red 500 (`#C74949`) — medium, matches small/medium nodes

## Favicon

The favicon is derived from `logo.svg`. At 32 × 32 px the three-level red scheme remains
distinguishable because the large anchor nodes (Red 900) contrast against the bond lines (Red 200).
