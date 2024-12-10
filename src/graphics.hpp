#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include "message.hpp"
#include "particle.hpp"

class Window {
    private:
        
        static const sf::Color brackgroundColor;
        static const sf::Color particleColor;
        static const double window_downscaler;
        static const int fps;

        sf::RenderWindow window;
        sf::Clock clock;
        int width;
        int height;
        double scale = 1;
        ParticleSet& particles;
        double radius;

        void draw(ParticleSet& particles);
        void display();

    public:
        Window(
            ParticleSet& particles,
            double radius = -1
        );

        void animate();
};


