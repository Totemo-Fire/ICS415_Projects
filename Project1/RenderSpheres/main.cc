#include <iostream>
#include "vector3.h"
#include <stdio.h>
#include "bmp.h"
#include "vector"

// Function Declarations

void canvasToViewPort(int, int, double*);
void intersectRaySphere(Vector3, Vector3, Vector3, double *);
Vector3 TraceRay(Vector3, Vector3, double, double);


// Set up the world (scene)

const int r = 1; // Sphere radius is 1
const Vector3 origin = Vector3();
Vector3 spheres[3];
Vector3 colors[4];

// Image and viewport properties

const int image_width = 256;
const int image_height = 256;
const double viewPortDistance = 1;
const double viewPortW = 1;
const double viewPortH = 1;

int main()
{
	// Generating spheres and colors
	spheres[0] = Vector3(0, -1, 3);
	spheres[1] = Vector3(2, 0, 4);
	spheres[2] = Vector3(-2, 0, 4);
	colors[0] = Vector3(255, 0, 0);
	colors[1] = Vector3(0, 0, 255);
	colors[2] = Vector3(0, 255, 0);
	colors[3] = Vector3(255, 255, 255);

	// Render

	std::vector<uint8_t> image(image_height * image_width * 3); // Create a buffer to hold the canvas data before writing it to the file
	for (int i = 0; i < image_width; i++)
	{
		for (int j = 0; j < image_height; j++)
		{
			double transformation[2]; // array to hold the result of transofrmation
			canvasToViewPort(i, j, transformation);
			Vector3 V = Vector3(transformation[0], transformation[1], viewPortDistance); // Calculate the V vector
			Vector3 color = TraceRay(origin, V + (-origin), 1, INFINITY); // Find the color of the pixel based on the ray and intersection with the spheres

			int index = (j * image_width + i) * 3; // Pixel location in the file

			image[index + 2] = (uint8_t) color.v[0];
			image[index + 1] = (uint8_t) color.v[1];
			image[index] = (uint8_t) color.v[2];
			
		}
	}
	writeBMP("image.bmp", image_width, image_height, image); // Finally write the canvas to the file
	return 0;
}

void canvasToViewPort(int i, int j, double* result)
{
	result[0] = (double) i * viewPortW / image_width -0.5; // -0.5 is the delta to transform the ceneter from the top left corner to the center
	result[1] = (double) j * viewPortH / image_height -0.5;
	return;
}

// res is the resultant t values
void intersectRaySphere(Vector3 O, Vector3 D, Vector3 sphere, double *res)
{
	Vector3 CO = O + (-sphere);
	double a = dot(D, D);
	double b = 2 * dot(CO, D);
	double c = dot(CO, CO) - r * r;
	double discriminant = b * b - 4 * a * c;
	if (discriminant < 0) // No solution
	{
		res[0] = INFINITY;
		res[1] = INFINITY;
		return;
	}
	double t1 = (-b + sqrt(discriminant)) / (2 * a);
	double t2 = (-b - sqrt(discriminant)) / (2 * a);
	res[0] = t1;
	res[1] = t2;
	return;
}

Vector3 TraceRay(Vector3 O, Vector3 D, double t_min, double t_max)
{
	double closest_t = INFINITY;
	int closest_sphere = 3; // 3 is used instead of null
	for (int i = 0; i < 3; i++)
	{
		double t[2]; // Store the t values
		intersectRaySphere(O, D, spheres[i], t);
		if ((t[0] < closest_t) && (t[0] >= t_min && t[0] < t_max)) // Ensure that the t value is at least a distance of 1 scene unit away from the camera, so that it can fit within the viewport
		{
			closest_t = t[0];
			closest_sphere = i;
		}
		if ((t[1] < closest_t) && (t[1] >= t_min && t[1] < t_max))
		{
			closest_t = t[1];
			closest_sphere = i;
		}
	}
	return colors[closest_sphere]; // Return the color based on the closest sphere determined 
}