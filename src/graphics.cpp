#include "graphics.hpp"


sf::Color const Window::brackgroundColor = sf::Color::Black;
sf::Color const Window::particleColor = sf::Color::White;
double const Window::window_downscaler = 1.4;
int const Window::fps = 10;




Window::Window(ParticleSet& particles, double radius) : particles(particles), radius(radius) {
    sf::VideoMode desktop = sf::VideoMode::getDesktopMode();
    width = desktop.width / window_downscaler;
    height = desktop.height / window_downscaler;
    window.create(sf::VideoMode(width, height), "N-Body Simulation");

    // let's put the window in the middle of the screen
    window.setPosition(sf::Vector2i(desktop.width / 2 - width / 2, desktop.height / 2 - height / 2));
}


void Window::display() {
    window.display();
}

void Window::draw(ParticleSet& particles) {
    window.clear(brackgroundColor);
    for (int i = 0; i < particles.size(); i++) {
        Particle p = particles.get(i);
        double screenX = p.position(0) * scale + width / 2;
        double screenY = - p.position(1) * scale + height / 2; // invert y axis

        sf::RectangleShape pixel(sf::Vector2f(1, 1)); // one pixel
        pixel.setFillColor(particleColor);
        pixel.setPosition(screenX, screenY); // position relative to top left corner 
        window.draw(pixel);
    }
}

void Window::animate() {

    std::vector<ParticleSet> sets = particles.split();
    radius = (radius==-1) ? sets[0].radius() * 1.1 : radius;
    radius = (radius == 0) ? 1 : radius;// avoid it to be 0

    scale = height / (2 * radius);
    int current_time = 0;
    double frameTime = (double) 1 / fps;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        current_time = current_time % sets.size();
        draw(sets[current_time]);
        display();
        current_time++;

        while (clock.getElapsedTime().asSeconds() < frameTime) {
            // Busy wait to maintain the frame rate
        }
        clock.restart();
        
    }
}