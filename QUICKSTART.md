# 🚀 Sagittarius A* Simulator - Quick Start

## Run in 3 Steps

### 1. Build (if not already done)
```bash
cd /home/dhh_elijah/black-hole-raytracer/build
make -j$(nproc) sagittarius-a
```

### 2. Run
```bash
./sagittarius-a
```

### 3. Explore!
- **W/S**: Camera distance
- **A/D**: Orbit around black hole
- **Mouse drag**: Free look
- **1-9**: Change black hole spin (try pressing 6, 7, 8, 9!)
- **I/K**: Adjust viewing angle

---

## 🎯 Cool Things to Try

### Demo 1: Watch the Disk Shrink
```
Press: 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9
Watch: Inner disk edge moves closer to horizon as spin increases!
```

### Demo 2: Doppler Asymmetry
```
1. Press '7' (high spin)
2. Press 'I' multiple times (edge-on view)
3. Observe: One side of disk is MUCH brighter (relativistic beaming!)
```

### Demo 3: Photon Ring
```
1. Press '9' (extremal spin)
2. Press 'W' to zoom in
3. Look for: Bright ring where photons orbit the black hole
```

### Demo 4: Fall Into the Black Hole
```
1. Press 'W' repeatedly to approach the event horizon
2. Watch: Disk wraps around you, starfield distorts
3. Cross r = 2M: Everything goes black (you're inside!)
```

---

## 📊 What the Display Shows

```
FPS: 87 | Spin: a*=0.9 | Incl: 17° | Camera: r=100
         ↑               ↑            ↑
         |               |            Distance from black hole
         |               Viewing angle (0° = face-on, 90° = edge-on)
         Current spin parameter (0.0 = Schwarzschild, 0.998 = extremal)
```

---

## ❓ Troubleshooting

### "ERROR: Could not find geodesic.comp shader!"
```bash
# Run from project root:
cd /home/dhh_elijah/black-hole-raytracer
./build/sagittarius-a
```

### Low FPS (<30)
- Your GPU may not support OpenGL 4.3+ compute shaders
- Check: `glxinfo | grep "OpenGL version"`
- Minimum: OpenGL 4.3 (released 2012)

### Black screen
- This is correct if you flew into the black hole (r < 2M)!
- Press 'S' to move back out

---

## 🎓 Learn More

- **Full documentation**: `SAGITTARIUS_A.md`
- **Technical details**: `REFACTOR_SUMMARY.md`
- **Controls**: Shown on startup in terminal

---

**Enjoy exploring the supermassive black hole at the center of our galaxy!** 🌌
