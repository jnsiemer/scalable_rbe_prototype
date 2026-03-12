// Copyright 2026 Jan Niklas Siemer
//
// This file is part of scalable_rbe.
//
// scalable_rbe is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::criterion_main;

mod rbe;

criterion_main! {rbe::benches}
