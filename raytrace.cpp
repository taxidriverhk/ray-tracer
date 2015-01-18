#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// Maximum depth of RGB calculation per pixel
const int MAX_DEPTH = 3;
// A very small number =/= 0
const float EPSILON = 0.0001;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.
struct Sphere
{
	string name;
	vec4 position; // Center
	vec4 scale; // Non-uniform scale
	vec4 color; // RGB values
	mat4 inverseTransformation; // Transformation matrix
	float Ka; // Constant for reflection of ambient light
	float Kd; // Constant for reflection of diffuse light
	float Ks; // Constant for reflection of specular light
	float Kr; // Constant for reflection
	float specularExponent;
	bool hollow; // Check if any part of the sphere cuts the near plane
				 // If so, then no illumination, shadow and reflection should be on the sphere
};

struct Light
{
	string name;
	vec4 position;
	vec4 intensity;
};

// Global variables
float g_near;
float g_left, g_right;
float g_top, g_bottom;

int g_width;
int g_height;

vector<vec4> g_colors;
vector<Sphere> g_spheres;
vector<Light> g_lights;

vec4 g_background_color;
vec4 g_ambient_intensity;

string g_output_filename;

// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    if(vs[0] == "NEAR")
		g_near = toFloat(vs[1]);
	else if(vs[0] == "LEFT")
		g_left = toFloat(vs[1]);
	else if(vs[0] == "RIGHT")
		g_right = toFloat(vs[1]);
	else if(vs[0] == "BOTTOM")
		g_bottom = toFloat(vs[1]);
	else if(vs[0] == "TOP")
		g_top = toFloat(vs[1]);
	else if(vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
	else if(vs[0] == "SPHERE")
	{
		Sphere sphere_buffer;
		
		// Parse data of the sphere
		sphere_buffer.name = vs[1];
		sphere_buffer.position = toVec4(vs[2], vs[3], vs[4]);
		sphere_buffer.scale = toVec4(vs[5], vs[6], vs[7]);
		sphere_buffer.color = toVec4(vs[8], vs[9], vs[10]);
		sphere_buffer.Ka = toFloat(vs[11]);
		sphere_buffer.Kd = toFloat(vs[12]);
		sphere_buffer.Ks = toFloat(vs[13]);
		sphere_buffer.Kr = toFloat(vs[14]);
		sphere_buffer.specularExponent = toFloat(vs[15]);
		
		// Check if the sphere cuts the near plane
		if(sphere_buffer.position.z + sphere_buffer.scale.z > -g_near
			&& sphere_buffer.position.z - sphere_buffer.scale.z < -g_near)
			sphere_buffer.hollow = true;
		else
			sphere_buffer.hollow = false;

		// Initialize the transformation matrix and invert it
		mat4 matrixBuffer;
		matrixBuffer *= Translate(sphere_buffer.position);
		matrixBuffer *= Scale(
			sphere_buffer.scale.x,
			sphere_buffer.scale.y,
			sphere_buffer.scale.z);
		InvertMatrix(matrixBuffer, sphere_buffer.inverseTransformation);

		// Add the sphere to the set of spheres 
		// (which should have up to 5 spheres)
		g_spheres.push_back(sphere_buffer);
	}
	else if(vs[0] == "LIGHT")
	{
		Light light_buffer;

		// Parse data of the light source
		light_buffer.name = vs[1];
		light_buffer.position = toVec4(vs[2], vs[3], vs[4]);
		light_buffer.intensity = toVec4(vs[5], vs[6], vs[7]);

		// Add the light to the set of lights
		// (which should have up to 5)
		g_lights.push_back(light_buffer);
	}
	else if(vs[0] == "BACK")
		g_background_color = toVec4(vs[1], vs[2], vs[3]);
	else if(vs[0] == "AMBIENT")
		g_ambient_intensity = toVec4(vs[1], vs[2], vs[3]);
	else if(vs[0] == "OUTPUT")
		g_output_filename = vs[1];
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
float getIntersection(const Ray& ray, const Sphere& sphere, bool shadowMode)
{
	// Untransform the ray function
	vec4 Sprime = sphere.inverseTransformation * ray.origin;
	vec4 cprime = sphere.inverseTransformation * ray.dir;

	// Solve for the quadratic equation
	float A = dot(cprime, cprime);
	float B = dot(Sprime, cprime);
	float C = dot(Sprime, Sprime) - 2;

	// If the ray intersects to a sphere, then the quadratic
	// equation above must have two solutions
	if((B*B - A*C) >= 0)
	{
		if(!sphere.hollow || shadowMode)
			return -B/A - sqrt(B*B - A*C)/A; // Take the smaller of the two solutions if the 
											 // sphere does not cut the near plane
		else
			return -B/A + sqrt(B*B - A*C)/A;
	}
	// Return 0 if the ray does not intersect the sphere or just tangent to it
	else
		return 0;
}

vec4 getNormal(const Ray& ray, const Sphere& sphere, float t)
{
	// Calculate the unit sphere function by untransforming it
	vec4 Sprime = sphere.inverseTransformation * ray.origin;
	vec4 cprime = sphere.inverseTransformation * ray.dir;

	// Get the hit point on the unit sphere
	vec4 unitSphereHitPoint = Sprime + cprime * t;
	// Get the normal of the hit point
	vec4 hitPointNormal = normalize(unitSphereHitPoint - vec4(0, 0, 0, 1));
	// Return the normal of the altered point
	vec4 alteredPointNormal = transpose(sphere.inverseTransformation) * hitPointNormal;
	alteredPointNormal.w = 0;
	return normalize(alteredPointNormal);
}

// -------------------------------------------------------------------
// Ray tracing

vec4 traceAuxiliary(const Ray& ray, int depth)
{
	if(depth == 0)
		return vec4(0, 0, 0, 1);

	// Initialize the resulting color
	vec4 resultColor(0, 0, 0, 1);

	// Find the closest intersection to any of the spheres
	float t_h = FLT_MAX;
	Sphere closestSphere;
	for(int i = 0; i < g_spheres.size(); i++)
	{
		// Try to find if the ray intersects to a certain unit sphere
		float t = getIntersection(ray, g_spheres[i], false);
		// Try to get the lowest t_h > 1, but reflected ray does not have to check such condition
		if(t < t_h && t > 1.0 && depth == MAX_DEPTH || (t > EPSILON && depth < MAX_DEPTH))
		{
			t_h = t;
			closestSphere = g_spheres[i];
		}
	}

	// No intersection is found, then color the pixel with background color
	if(t_h == FLT_MAX)
	{
		if(depth == MAX_DEPTH)
			return g_background_color;
		// For reflected ray
		else if(depth < MAX_DEPTH)
			return vec4(0, 0, 0, 1);
	}

	// Get the hit point to the transformed sphere
	vec4 hitPoint = ray.origin + ray.dir * t_h;
	vec4 normal = getNormal(ray, closestSphere, t_h);

	// Calculate the resulting intensity
	for(int j = 0; j < g_lights.size(); j++)
	{
		// Calculate the light, view and reflection vectors, which are unit vectors
		vec4 light = normalize(g_lights[j].position - hitPoint);
		vec4 view = normalize(ray.origin - hitPoint);
		vec4 reflection = normalize(2*normal*dot(normal, light) - light);
		light.w = 0; view.w = 0; reflection.w = 0;

		// Check if the light source shadows any of the spheres
		bool useLight = true;
		// Set up a ray from hit point to light source
		Ray hitPointRay;
		hitPointRay.origin = hitPoint;
		hitPointRay.dir = g_lights[j].position - hitPoint;
		for(int k = 0; k < g_spheres.size(); k++)
		{
			// The function already has S' and c' calculated
			float t_min = getIntersection(hitPointRay, g_spheres[k], true);
			// If there are any intersections between 0.0001 and 1, 
			// then break the loop and don't use light
			if(t_min > EPSILON && t_min < 1)
			{
				useLight = false;
				break;
			}
		}

		// Add the intensity only if the angle between normal and light <= 90 degrees
		// , because if their angle > 90, then the eye can not be able to see the light
		// Also, if no sphere intersects the light source, then don't let light
		// contribute to the specular and diffuse
		if(dot(normal, light) > 0 && useLight)
		{
			// Phong illumination model
			vec4 resultIntensity;
			resultIntensity += closestSphere.Kd * g_lights[j].intensity * dot(normal, light) * closestSphere.color; // Diffuse
			// Ignore specular part if sphere is hollow
			if(!closestSphere.hollow)
				resultIntensity	+= closestSphere.Ks * g_lights[j].intensity * pow(dot(reflection, view), closestSphere.specularExponent); // Specular
			resultColor += resultIntensity;
		}
	}


	Ray reflectedRay;
	// If the angle between the normal and ray is less than 90 degrees
	// , then the eye would not be able to see the reflection
	if(dot(normal, ray.dir) < 0)
	{
		reflectedRay.origin = hitPoint;
		reflectedRay.dir = normalize(-2*dot(normal, ray.dir)*normal + ray.dir);
	}
	else
	{
		reflectedRay.origin = vec4(0, 0, 0, 1);
		reflectedRay.dir = vec4(0, 0, 0, 0);
	}

	resultColor += closestSphere.Ka * g_ambient_intensity * closestSphere.color 
		+ closestSphere.Kr * traceAuxiliary(reflectedRay, depth-1);
	return resultColor;
}

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
	return traceAuxiliary(ray, MAX_DEPTH);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    return vec4(
		// Get the actual dimension first, then get actual position by multiplying the percentage
		g_left + ((float)ix)/g_width * (g_right - g_left), 
		g_bottom + ((float)iy)/g_height * (g_top - g_bottom), 
		-g_near, 
		0.0);
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);

	// Clamp values if out of range
    if(color.x > 1.0) color.x = 1.0;
	if(color.y > 1.0) color.y = 1.0;
	if(color.z > 1.0) color.z = 1.0;

	setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
			{
				if(g_colors[y*g_width+x].x > 1.0) g_colors[y*g_width+x].x = 1.0;
				if(g_colors[y*g_width+x].y > 1.0) g_colors[y*g_width+x].y = 1.0;
				if(g_colors[y*g_width+x].z > 1.0) g_colors[y*g_width+x].z = 1.0;
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
			}

    // TODO: change file name based on input file name.
	char* temp_filename = &g_output_filename[0];
    savePPM(g_width, g_height, temp_filename, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    //if (argc < 2)
    //{
    //    cout << "Usage: template-rt <input_file.txt>" << endl;
    //    exit(1);
    //}
	loadFile("testIllum.txt");
	render();
	saveFile();
	return 0;
}

