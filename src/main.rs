pub mod cannon;
pub mod canvas;
pub mod color;
pub mod ppm_file;
pub mod tuple;
pub mod matrix;
pub mod transformation;

#[macro_use]
extern crate approx;

use ppm_file::PpmFile;
use tuple::Normalize;

use crate::cannon::{environment, projectile, tick};
use crate::canvas::canvas;
use crate::color::Colors;
use crate::tuple::{point, vector};

fn main() {
    let start = point(0, 1, 0);
    let velocity = vector(1.0, 1.8, 0.0).normalize() * 11.25;
    let mut bullet = projectile(start.into(), velocity);
    let gravity = vector(0.0, -0.1, 0.0);
    let wind = vector(-0.01, 0.0, 0.0);
    let env = environment(gravity, wind);

    let mut canv = canvas(900, 500);
    let red = Colors::red();

    while bullet.on_surface() {
        bullet = tick(&env, &bullet);
        let position = &bullet.position;
        let x = position.x as usize;
        let y = position.y as usize;

        canv.write_pixel(x, 500 - y, red);
    }

    let file = PpmFile::from_canvas(&canv);
    file.save("output/canvas.ppm").unwrap_or_default();
}
