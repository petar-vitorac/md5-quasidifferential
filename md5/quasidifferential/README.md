# MD5 - Quasidifferential

- `requirements.txt` contains the list of Python dependencies. They can be installed with `pip install -r requirements.txt`.
- `qd-trails.py` computes linearly independent quasidifferential trails in MD5 for the given characteristic with intermediate differences. To only find deterministic trails, first uncomment `assert_masks_0()` to find the weight when the masks are 0. Then, change `max_target_weight` so that it matches the weight of the deterministic trails (same weight as when masks are 0).
