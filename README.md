# Classical sine-Gordon quench

This repository now contains a static browser demo for a closed-loop sine-Gordon simulation.

## Files

- `index.html` — page layout and controls
- `style.css` — styling
- `sim.js` — periodic sine-Gordon solver and canvas rendering

## Model

The demo evolves the periodic sine-Gordon equation

`phi_tt - phi_xx + sin(phi) = 0`

with initial condition

`phi(x,0) = 0`

`phi_t(x,0) = B exp(-(x/sigma)^2)`

on a ring of circumference `L`.

## Running locally

Open `index.html` in a browser, or serve the folder with a static file server.

## GitHub Pages

To publish this as a website, enable **Settings → Pages** and choose the `main` branch as the source.
