#include <iostream>
#include "vector3.h"
#include "bmp.h"

// Function Declarations

void canvasToViewPort(int, int, double *);
void intersectRaySphere(Vector3, Vector3, int, double *);
Vector3 TraceRay(Vector3, Vector3, double, double);
double ComputeLighting(Vector3, Vector3, Vector3, double);
Vector3 applyColorRangeCeiling(Vector3);

// Set up the world (scene)

const int radii[4] = {1, 1, 1, 5000}; // Sphere radius is 1
const Vector3 origin = Vector3();
Vector3 spheres[4];
Vector3 colors[5];
double ambientLight = 0.2; // One ambient light for the whole scene
double pointLights[1];	   // Point lights
Vector3 pointLightsPositions[1];
double directionalLights[1];
Vector3 directionalLightsDirections[1];
double specular[4];

// Image and viewport properties

const int image_width = 512;
const int image_height = 512;
const double viewPortDistance = 1;
const double viewPortW = 1;
const double viewPortH = 1;

int main()
{
	// Generating spheres, light and colors
	spheres[0] = Vector3(0, -1, 3);
	spheres[1] = Vector3(2, 0, 4);
	spheres[2] = Vector3(-2, 0, 4);
	spheres[3] = Vector3(0, -5001, 0);
	colors[0] = Vector3(255, 0, 0);
	colors[1] = Vector3(0, 0, 255);
	colors[2] = Vector3(0, 255, 0);
	colors[3] = Vector3(255, 255, 0);
	colors[4] = Vector3(255, 255, 255);
	pointLights[0] = 0.6;
	pointLightsPositions[0] = Vector3(2, 1, 0);
	directionalLights[0] = 0.2;
	directionalLightsDirections[0] = Vector3(1, 4, 4);
	specular[0] = 500;
	specular[1] = 500;
	specular[2] = 10;
	specular[3] = 5000;

	bool DEBUGGING = 0;
	if (DEBUGGING) {
		int x = 254;
		int y = 0;
		double cv[2];
		canvasToViewPort(x, y, cv);
		printf("x: %d, y: %d, cvi: %lf, cvj: %lf\n", x, y, cv[0], cv[1]);
		Vector3 V = Vector3(cv[0], cv[1], viewPortDistance);
		Vector3 color = TraceRay(origin, V +(-origin), 1, INFINITY);
		printf("V: %lf, %lf, %lf\nTraceRay color: %lf, %lf, %lf\n", V.v[0], V.v[1], V.v[2], color.v[0], color.v[1], color.v[2]);
		return 0;
	}


	// Render

	std::vector<uint8_t> image(image_height * image_width * 3); // Create a buffer to hold the canvas data before writing it to the file
	for (int i = 0; i < image_width; i++)
	{
		for (int j = 0; j < image_height; j++)
		{
			double transformation[2]; // array to hold the result of transofrmation
			canvasToViewPort(i, j, transformation);
			Vector3 V = Vector3(transformation[0], transformation[1], viewPortDistance); // Calculate the V vector
			Vector3 color = TraceRay(origin, V + (-origin), 1, INFINITY);				 // Find the color of the pixel based on the ray and intersection with the spheres and light
			int index = (j * image_width + i) * 3; // Pixel location in the file
			image[index + 2] = (uint8_t)color.v[0];
			image[index + 1] = (uint8_t)color.v[1];
			image[index] = (uint8_t)color.v[2];
		}
	}
	writeBMP("image.bmp", image_width, image_height, image); // Finally write the canvas to the file
	return 0;
}

void canvasToViewPort(int i, int j, double *result)
{
	result[0] = (double)i * viewPortW / image_width - 0.5; // -0.5 is the delta to transform the canvas from the top left corner to the center
	result[1] = (double)j * viewPortH / image_height - 0.5;
	return;
}

// res is the resultant t values
void intersectRaySphere(Vector3 O, Vector3 D, int sphereID, double *res)
{
	Vector3 CO = O + (-spheres[sphereID]);
	double a = dot(D, D);
	double b = 2 * dot(CO, D);
	double c = dot(CO, CO) - radii[sphereID] * radii[sphereID];
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
	int closest_sphere = -1; // -1 is used instead of null
	for (int i = 0; i < 4; i++) // Loop through all spheres
	{
		double t[2]; // Store the t values
		intersectRaySphere(O, D, i, t);
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
	if(closest_sphere == -1) {
		return colors[4];
	}
	Vector3 P = O + D * closest_t;
	Vector3 N = P + -spheres[closest_sphere];
	N = N / (N.magnitude());
	Vector3 result = colors[closest_sphere] * ComputeLighting(P, N, -D, specular[closest_sphere]); // Return the color based on the closest sphere determined and lighting
	return applyColorRangeCeiling(result); // Ensure that color values stay at or below 255;
}

double ComputeLighting(Vector3 P, Vector3 N, Vector3 V, double s)
{	
	double intensity = ambientLight;
	Vector3 L;
	double n_dot_l;
	for (int i = 0; i < 1; i++)
	{ // Loop through all point lights
		L = pointLightsPositions[i] + (-P);
		n_dot_l = dot(N, L);

		// Diffusion
		if (n_dot_l > 0)
		{
			intensity += pointLights[i] * n_dot_l / (N.magnitude() * L.magnitude());
		}

		// Specular
		if (s >= 0)
		{
			Vector3 R = N * 2 * n_dot_l + (-L);
			double r_dot_v = dot(R, V);
			if (r_dot_v > 0)
			{
				intensity += pointLights[i] * pow(r_dot_v/(R.magnitude() * V.magnitude()), s);
			}
		}
	}

	for (int i = 0; i < 1; i++)
	{ // Loop through all directional lights
		L = directionalLightsDirections[i];
		
		// Diffusion
		n_dot_l = dot(N, L);
		if (n_dot_l > 0)
		{
			intensity += directionalLights[i] * n_dot_l / (N.magnitude() * L.magnitude());
		}

		// Specular
		if (s >= 0)
		{
			Vector3 R = N * 2 * n_dot_l + (-L);
			double r_dot_v = dot(R, V);
			if (r_dot_v > 0)
			{
				intensity += directionalLights[i] * pow(r_dot_v/(R.magnitude() * V.magnitude()), s);
			}
		}
	}
	return intensity;
}

Vector3 applyColorRangeCeiling(Vector3 u)
{
	return Vector3(fmin(255, u.v[0]), fmin(255, u.v[1]), fmin(255, u.v[2]));
}