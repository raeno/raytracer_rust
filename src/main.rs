pub mod cannon;
pub mod canvas;
pub mod color;
pub mod matrix;
pub mod ppm_file;
pub mod transformation;
pub mod tuple;

#[macro_use]
#[cfg(test)]
extern crate approx;

use ppm_file::PpmFile;
use tuple::Normalize;

use crate::cannon::{environment, projectile, tick};
use crate::canvas::{Canvas, canvas};
use crate::color::Colors;
use crate::tuple::{point, vector};

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "Raytracer Challenge")]
#[command(version = "0.4")]
#[command(about = "Solutions from raytracer challenge book", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Cannon,
    Clock,
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Cannon => {
            cannon_game();
        },
        Commands::Clock => {
            clock_game()
        }
    }
}

fn cannon_game() {
    save_to_file(&fire_cannon(), "output/cannon.ppm");
}

fn clock_game() {
    save_to_file(
        &draw_clock(),
        "output/clock.ppm"
    )
}

fn fire_cannon() -> Canvas {
    let start = point(0, 1, 0);
    let velocity = vector(1.0, 1.8, 0.0).normalize() * 11.25;
    let mut bullet = projectile(start.into(), velocity);
    let gravity = vector(0.0, -0.1, 0.0);
    let wind = vector(-0.01, 0.0, 0.0);
    let env = environment(gravity, wind);

    let mut canvas = canvas(900, 500);

    while bullet.on_surface() {
        bullet = tick(&env, &bullet);
        let position = &bullet.position;
        let x = position.x as usize;
        let y = position.y as usize;

        canvas.write_pixel(x, 500 - y, Colors::red());
    }
    canvas
}

fn draw_clock() -> Canvas {
    let canvas = canvas(900, 500);
    canvas
}

fn save_to_file(canvas: &Canvas, path: &str) {
    let file = PpmFile::from_canvas(&canvas);
    file.save(path).unwrap_or_default();
}
